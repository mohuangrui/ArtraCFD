/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>                        *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "immersed_boundary.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <float.h> /* size of floating point values */
#include "computational_geometry.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    LAYER1 = 6,
    LAYER2 = 24,
    LAYER3 = 30,
    R = 2,
} IBM;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void InitializeGeometryDomain(Space *);
static void IdentifyGeometryNode(Space *);
static void IdentifyInterfacialNode(Space *, const Model *);
static int InterfacialState(const int, const int, const int, const int, const int,
        const int [restrict][DIMS], const Node *const , const Partition *);
static int GhostState(const int, const int, const int, const int, const int,
        const int [restrict][DIMS], const Node *const , const Partition *);
static void ApplyWeighting(const Real [restrict], const Real, Real, 
        Real [restrict], Real [restrict]);
static Real InverseDistanceWeighting(const int, const int [restrict], const Real [restrict],
        const int, const int, const int, const Partition *, const Node *const , const Model *, Real [restrict]);
static void FlowReconstruction(const int, const int [restrict], const Real [restrict], const int,
        const int, const int, const Polyhedron *, const Partition *, const Node *const , const Model *,
        const Real [restrict], const Real [restrict], Real [restrict], Real [restrict]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * It's necessary to differentiate exterior nodes and interior nodes to avoid
 * incorrect mark of interfacial nodes near global domain boundaries.
 *
 * The identification process should proceed step by step to correctly handle
 * all the relationships and avoid interference between each step,
 * interference may happen if identification processes are crunched together.
 *
 * Moreover, whenever identifying a node, link the node to geometry, which 
 * can be used to access the geometric data for computing required data. 
 * The rational is that don't store every information for each node, but
 * only store necessary information. When need it, access and calculate it.
 */
void ComputeGeometryDomain(Space *space, const Model *model)
{
    InitializeGeometryDomain(space);
    IdentifyGeometryNode(space);
    IdentifyInterfacialNode(space, model);
    return;
}
static void InitializeGeometryDomain(Space *space)
{
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                node[idx].geoID = 0;
                node[idx].layerID = 0;
                node[idx].ghostID = 0;
            }
        }
    }
    return;
}
/*
 * When locate geometry nodes, there are two approaches available. One is search
 * over each node and verify each node regarding to all the geometries; another
 * is search each geometry and find all the nodes inside current geometry.
 * The second method is adopted here for performance reason.
 *
 * Points either in or on geometry should be classified into the corresponding
 * geometry, which is to best utilize the convergence property of immersed
 * boundary treatment.
 *
 * It is efficient to only test points that are inside the bounding box or
 * sphere of a large polyhedron. Be cautious with the validity of any calculated
 * index. It's extremely necessary to adjust the index into the valid region or 
 * check the validity of the index to avoid index exceed array bound limits.
 *
 * For a large data set, an extra data preprocessing with spatial subdivision
 * is required to ensure minimal test in addition to the bounding container
 * method. Spatial subdivision is to provide internal resolution for the
 * polyhedron for fast inclusion determination.
 */
static void IdentifyGeometryNode(Space *space)
{
    const Partition *restrict part = &(space->part);
    Geometry *geo = &(space->geo);
    Node *const node = space->node;
    Polyhedron *poly = NULL;
    const IntVec nMin = {part->ns[PIN][X][MIN], part->ns[PIN][Y][MIN], part->ns[PIN][Z][MIN]};
    const IntVec nMax = {part->ns[PIN][X][MAX], part->ns[PIN][Y][MAX], part->ns[PIN][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const int ng = part->ng;
    int faceID = 0;
    int idx = 0; /* linear array index math variable */
    int box[DIMS][LIMIT] = {{0}}; /* bounding box in node space */
    RealVec p = {0.0}; /* node point */
    for (int n = 0; n < geo->totalN; ++n) {
        poly = geo->poly + n;
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ValidNodeSpace(NodeSpace(poly->box[s][MIN], sMin[s], dd[s], ng), nMin[s], nMax[s]);
            box[s][MAX] = ValidNodeSpace(NodeSpace(poly->box[s][MAX], sMin[s], dd[s], ng), nMin[s], nMax[s]) + 1;
        }
        /* find nodes in geometry, then flag and link to geometry. */
        for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
            for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    if (0 != node[idx].geoID) { /* already identified by others */
                        continue;
                    }
                    p[X] = PointSpace(i, sMin[X], d[X], ng);
                    p[Y] = PointSpace(j, sMin[Y], d[Y], ng);
                    p[Z] = PointSpace(k, sMin[Z], d[Z], ng);
                    if (0 == poly->faceN) { /* analytical sphere */
                        if (poly->r * poly->r >= Dist2(poly->O, p)) {
                            node[idx].geoID = n + 1;
                            node[idx].faceID = 0;
                        }
                    } else { /* triangulated polyhedron */
                        if (PointInPolyhedron(p, poly, &faceID)) {
                            node[idx].geoID = n + 1;
                            node[idx].faceID = faceID;
                        }
                    }
                }
            }
        }
    }
    return;
}
/*
 * Convective terms only have first order derivatives, discretization
 * stencils of each computational nodes are cross types.
 * Diffusive terms have second order mixed derivatives, discretization
 * stencils of each computational nodes are plane squares with conner
 * nodes. Therefore, for interfacial computational nodes, neighbouring
 * nodes at corner directions require to be considered and is dependent
 * on the discretization of diffusive fluxes. An interfacial node is a
 * node that has heterogeneous neighbouring nodes on the searching path.
 * A ghost node is a node that is outside the normal computational domain
 * but locates on the numerical boundary. Hence, ghost nodes form a subset
 * of the interfacial nodes.
 */
/*
 * When dealing with moving geometries, interfacial nodes may change their
 * corresponding geometry to others and keep as interfacial nodes in the 
 * new geometry. When a no ghost approach is adopted, this change will not
 * bring any problems since for each geometry, the newly joined nodes always
 * become to boundary nodes and their values will be constructed by boundary
 * treatment. However, for a ghost approach, the newly jointed nodes directly
 * become to normal computational nodes in the new geometry, therefore, a 
 * reconstruction is required to correctly deal with these newly jointed nodes.
 * To detect these nodes, two states, geometry identifier and face identifier 
 * are used together to flag nodes in geometry. After mesh regeneration, the
 * geometry identifier is updated and the face identifier maintains old value.
 * By comparing the updated value with the old value, we then are able to detect
 * those nodes previously are not in the current geometry and now newly become 
 * normal nodes in the current geometry. The physical values of the detected
 * nodes should be reconstructed from normal nodes with a NONE face identifier.
 * After that, the face identifier shall also be reset to NONE.
 * Note: if to apply this mechanism to all the regions, a new field identifier
 * that copy the geometry identifier shall be used rather than simply reusing
 * the face identifier to play this dual role.
 * Note: for two area touching objects separating instantly, fresh normal
 * nodes are generated without any valid neighbours, which situation requires 
 * further treatment, which is unsolved currently.
 */
static void IdentifyInterfacialNode(Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    const int end = LAYER2 + (part->ng - 1) * (LAYER3 - LAYER2); /* max search layer is ng + 1 */
    IntVec n = {0};
    RealVec p = {0.0};
    Real Uo[DIMUo] = {0.0};
    Real weightSum = 0.0;
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                if ((NONE != node[idx].faceID) && (0 == node[idx].geoID)) {
                    /* a newly joined node */
                    n[X] = i;
                    n[Y] = j;
                    n[Z] = k;
                    p[X] = PointSpace(i, part->domain[X][MIN], part->d[X], part->ng);
                    p[Y] = PointSpace(j, part->domain[Y][MIN], part->d[Y], part->ng);
                    p[Z] = PointSpace(k, part->domain[Z][MIN], part->d[Z], part->ng);
                    weightSum = InverseDistanceWeighting(TO, n, p, R, NONE, 0, part, node, model, Uo);
                    /* Normalize the weighted values */
                    Normalize(DIMUo, weightSum, Uo);
                    Uo[0] = Uo[4] / (Uo[5] * model->gasR); /* compute density */
                    ConservativeByPrimitive(model->gamma, Uo, node[idx].U[TO]);
                    node[idx].faceID = NONE; /* reset after correct reconstruction */
                }
                if (0 == node[idx].geoID) { /* skip interfacial nodes for main domain */
                    continue;
                }
                /* search neighbours to determine interfacial state */
                node[idx].layerID = InterfacialState(k, j, i, node[idx].geoID, end, part->path, node, part);
                if ((0 != node[idx].layerID) && (0 != node[idx].geoID)) { /* an interfacial node may be a ghost node */
                    /* search neighbours to determine ghost state */
                    node[idx].ghostID = GhostState(k, j, i, 0, end, part->path, node, part);
                }
            }
        }
    }
    return;
}
static int InterfacialState(const int k, const int j, const int i, const int geoID, const int end,
        const int path[restrict][DIMS], const Node *const node, const Partition *part)
{
    /*
     * Search around the specified node to check whether current node
     * is an interfacial node, and return the interfacial state.
     */
    int idx = 0; /* linear array index math variable */
    for (int n = 0; n < end; ++n) {
        idx = IndexNode(k + path[n][Z], j + path[n][Y], i + path[n][X], part->n[Y], part->n[X]);
        if (NONE == node[idx].geoID) { /* an exterior node is not valid */
            continue;
        }
        if (geoID != node[idx].geoID) { /* a heterogeneous node on the path */
            if (LAYER1 > n) {
                return 1;
            } else {
                if (LAYER2 > n) {
                    return 2;
                } else {
                    return 3;
                }
            }
        }
    }
    return 0;
}
static int GhostState(const int k, const int j, const int i, const int geoID, const int end,
        const int path[restrict][DIMS], const Node *const node, const Partition *part)
{
    /*
     * Search around the specified node to check whether current node
     * is on numerical boundary, and return the state.
     */
    int idx = 0; /* linear array index math variable */
    for (int n = 0; n < end; ++n) {
        idx = IndexNode(k + path[n][Z], j + path[n][Y], i + path[n][X], part->n[Y], part->n[X]);
        if (NONE == node[idx].geoID) { /* an exterior node is not valid */
            continue;
        }
        if (geoID == node[idx].geoID) { /* a normal computational node on the path */
            if (LAYER1 > n) {
                return 1;
            } else {
                if (LAYER2 > n) {
                    return 2;
                } else {
                    return 3;
                }
            }
        }
    }
    return 0;
}
void ImmersedBoundaryTreatment(const int tn, Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    const IntVec nMin = {part->ns[PIN][X][MIN], part->ns[PIN][Y][MIN], part->ns[PIN][Z][MIN]};
    const IntVec nMax = {part->ns[PIN][X][MAX], part->ns[PIN][Y][MAX], part->ns[PIN][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const int ng = part->ng;
    IntVec nI = {0}; /* image node */
    IntVec nG = {0}; /* ghost node */
    RealVec pG = {0.0}; /* ghost point */
    RealVec pO = {0.0}; /* boundary point */
    RealVec pI = {0.0}; /* image point */
    RealVec N = {0.0}; /* normal */
    Real UoG[DIMUo] = {0.0};
    Real UoO[DIMUo] = {0.0};
    Real UoI[DIMUo] = {0.0};
    Real weightSum = 0.0;
    int box[DIMS][LIMIT] = {{0}}; /* bounding box in node space */
    for (int n = 0; n < geo->totalN; ++n) {
        poly = geo->poly + n;
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ValidNodeSpace(NodeSpace(poly->box[s][MIN], sMin[s], dd[s], ng), nMin[s], nMax[s]);
            box[s][MAX] = ValidNodeSpace(NodeSpace(poly->box[s][MAX], sMin[s], dd[s], ng), nMin[s], nMax[s]) + 1;
        }
    for (int r = 1; r < ng + 2; ++r) { /* layer by layer treatment */
        for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
            for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    if ((r != node[idx].ghostID) || (n + 1 != node[idx].geoID)) {
                        continue;
                    }
                    pG[X] = PointSpace(i, sMin[X], d[X], ng);
                    pG[Y] = PointSpace(j, sMin[Y], d[Y], ng);
                    pG[Z] = PointSpace(k, sMin[Z], d[Z], ng);
                    if (model->layers >= r) { /* immersed boundary treatment */
                        ComputeGeometricData(node[idx].faceID, poly, pG, pO, pI, N);
                        nI[X] = NodeSpace(pI[X], sMin[X], dd[X], ng);
                        nI[Y] = NodeSpace(pI[Y], sMin[Y], dd[Y], ng);
                        nI[Z] = NodeSpace(pI[Z], sMin[Z], dd[Z], ng);
                        /*
                         * When extremely strong discontinuities exist in the
                         * domain of dependence of inverse distance weighting,
                         * the original weighting approach may produce spurious
                         * discontinuities among the weighted node and its
                         * neighboring nodes. The best way to solve this issue
                         * is to apply WENO's idea to avoid discontinuous
                         * stencils and to only use smooth stencils. However,
                         * this will make the algorithm much more complex.
                         */
                        FlowReconstruction(tn, nI, pI, R, NONE, 0, poly, part, node, model, pO, N, UoO, UoI);
                        MethodOfImage(UoI, UoO, UoG);
                    } else { /* inverse distance weighting */
                        nG[X] = i;
                        nG[Y] = j;
                        nG[Z] = k;
                        weightSum = InverseDistanceWeighting(tn, nG, pG, 1, r - 1, n + 1, part, node, model, UoG);
                        /* Normalize the weighted values */
                        Normalize(DIMUo, weightSum, UoG);
                    }
                    UoG[0] = UoG[4] / (UoG[5] * model->gasR); /* compute density */
                    ConservativeByPrimitive(model->gamma, UoG, node[idx].U[tn]);
                }
            }
        }
    }
    }
    return;
}
void MethodOfImage(const Real UoI[restrict], const Real UoO[restrict], Real UoG[restrict])
{
    /* 
     * Apply the method of image.
     *  -- reflecting vectors over wall for both slip and noslip, stationary and
     *     moving conditions is unified by linear interpolation.
     *  -- scalars are symmetrically reflected between image and ghost.
     *  -- other scalars are determined by equation of state.
     */
    UoG[1] = 2.0 * UoO[1] - UoI[1];
    UoG[2] = 2.0 * UoO[2] - UoI[2];
    UoG[3] = 2.0 * UoO[3] - UoI[3];
    UoG[4] = UoI[4];
    UoG[5] = UoI[5];
    return;
}
void ComputeGeometricData(const int faceID, const Polyhedron *poly, const Real pG[restrict],
        Real pO[restrict], Real pI[restrict], Real N[restrict])
{
    if (0 == poly->faceN) { /* analytical sphere */
        Real dist = 0.0;
        N[X] = pG[X] - poly->O[X];
        N[Y] = pG[Y] - poly->O[Y];
        N[Z] = pG[Z] - poly->O[Z];
        dist = Norm(N);
        Normalize(DIMS, dist, N);
        dist = poly->r - dist;
        pO[X] = pG[X] + dist * N[X];
        pO[Y] = pG[Y] + dist * N[Y];
        pO[Z] = pG[Z] + dist * N[Z];
    } else { /* triangulated polyhedron */
        ComputeIntersection(pG, faceID, poly, pO, N);
    }
    pI[X] = pO[X] + pO[X] - pG[X];
    pI[Y] = pO[Y] + pO[Y] - pG[Y];
    pI[Z] = pO[Z] + pO[Z] - pG[Z];
    return;
}
static void FlowReconstruction(const int tn, const int n[restrict], const Real p[restrict], const int h,
        const int type, const int geoID, const Polyhedron *poly, const Partition *part, const Node *const node, 
        const Model *model, const Real pO[restrict], const Real N[restrict], Real UoO[restrict], Real UoI[restrict])
{
    Real weightSum = 0.0; /* store the sum of weights */
    /* pre-estimate step */
    weightSum = InverseDistanceWeighting(tn, n, p, h, type, geoID, part, node, model, UoI);
    /* physical boundary condition enforcement step */
    RealVec Vs = {0.0}; /* general motion of boundary point */
    /* Vs = Vcentroid + W x r */
    const RealVec r = {pO[X] - poly->O[X], pO[Y] - poly->O[Y], pO[Z] - poly->O[Z]};
    Cross(poly->W, r, Vs); /* relative motion in translating coordinate system */
    Vs[X] = poly->V[X] + Vs[X];
    Vs[Y] = poly->V[Y] + Vs[Y];
    Vs[Z] = poly->V[Z] + Vs[Z];
    if (0 < poly->cf) { /* noslip wall */
        UoO[1] = Vs[X];
        UoO[2] = Vs[Y];
        UoO[3] = Vs[Z];
    } else { /* slip wall */
        RealVec Ta = {0.0}; /* tangential vector */
        RealVec Tb = {0.0}; /* tangential vector */
        Real RHS[DIMS] = {0.0}; /* right hand side vector */
        OrthogonalSpace(N, Ta, Tb);
        RHS[X] = Dot(Vs, N);
        RHS[Y] = Dot(UoI + 1, Ta) / weightSum;
        RHS[Z] = Dot(UoI + 1, Tb) / weightSum;
        UoO[1] = N[X] * RHS[X] + Ta[X] * RHS[Y] + Tb[X] * RHS[Z];
        UoO[2] = N[Y] * RHS[X] + Ta[Y] * RHS[Y] + Tb[Y] * RHS[Z];
        UoO[3] = N[Z] * RHS[X] + Ta[Z] * RHS[Y] + Tb[Z] * RHS[Z];
    }
    UoO[4] = UoI[4] / weightSum; /* zero gradient pressure */
    if (0.0 > poly->T) { /* adiabatic, dT/dn = 0 */
        UoO[5] = UoI[5] / weightSum;
    } else { /* otherwise, use specified constant wall temperature, T = Tw */
        UoO[5] = poly->T;
    }
    UoO[0] = UoO[4] / (UoO[5] * model->gasR); /* compute density */
    /* correction step by adding the boundary point as a stencil */
    ApplyWeighting(UoO, part->tinyL, Dist2(p, pO), &weightSum, UoI);
    /* Normalize the weighted values */
    Normalize(DIMUo, weightSum, UoI);
    return;
}
static Real InverseDistanceWeighting(const int tn, const int n[restrict], const Real p[restrict], 
        const int h, const int type, const int geoID, const Partition *part, 
        const Node *const node, const Model *model, Real Uo[restrict])
{
    int idx = 0; /* linear array index math variable */
    const int idxMax = part->n[X] * part->n[Y] * part->n[Z];
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const int ng = part->ng;
    int tally = 0; /* stencil count and zero stencil detector */
    Real Uoh[DIMUo] = {0.0}; /* primitive at neighbouring node */
    RealVec ph = {0.0}; /* neighbouring point */
    Real weightSum = 0.0;
    for (int dim = 0; dim < DIMUo; ++dim) {
        Uo[dim] = 0.0; /* reset */
    }
    /* 
     * Search nodes with required "type" in the domain specified by the center node
     * "n" and range "h" as interpolation stencils for the interpolated point "p".
     */
    for (int kh = -h; kh <= h; ++kh) {
        for (int jh = -h; jh <= h; ++jh) {
            for (int ih = -h; ih <= h; ++ih) {
                idx = IndexNode(n[Z] + kh, n[Y] + jh, n[X] + ih, part->n[Y], part->n[X]);
                if ((0 > idx) || (idxMax <= idx)) {
                    continue; /* illegal index */
                }
                if (geoID != node[idx].geoID) {
                    continue;
                }
                if (0 == geoID) { /* require normal node type */
                    if (type != node[idx].faceID) {
                        continue;
                    }
                } else { /* require specified ghost node type */
                    if (type != node[idx].ghostID) { /* not a ghost node with current type */
                        continue;
                    }
                }
                ++tally;
                ph[X] = PointSpace(n[X] + ih, sMin[X], d[X], ng);
                ph[Y] = PointSpace(n[Y] + jh, sMin[Y], d[Y], ng);
                ph[Z] = PointSpace(n[Z] + kh, sMin[Z], d[Z], ng);
                PrimitiveByConservative(model->gamma, model->gasR, node[idx].U[tn], Uoh);
                /* use distance square to avoid expensive sqrt */
                ApplyWeighting(Uoh, part->tinyL, Dist2(p, ph), &weightSum, Uo);
            }
        }
    }
    if (0 == tally) {
        fprintf(stderr, "k=%d, j=%d, i=%d\n", n[Z], n[Y], n[X]);
        FatalError("zero stencil encountered...");
    }
    return weightSum;
}
static void ApplyWeighting(const Real Uoh[restrict], const Real tiny, Real weight, 
        Real weightSum[restrict], Real Uo[restrict])
{
    if (tiny > weight) { /* avoid overflow of too small weight */
        weight = tiny;
    }
    weight = 1.0 / weight; /* compute weight */
    for (int n = 0; n < DIMUo; ++n) {
        Uo[n] = Uo[n] + Uoh[n] * weight;
    }
    *weightSum = *weightSum + weight; /* accumulate normalizer */
    return;
}
/* a good practice: end file with a newline */

