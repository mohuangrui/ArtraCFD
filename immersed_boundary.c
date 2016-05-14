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
    LAYER1 = 18,
    LAYER2 = 24,
    LAYER3 = 30,
} PATH;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void InitializeGeometryDomain(Space *);
static void IdentifyGeometryNode(Space *);
static void IdentifyInterfacialNode(Space *);
static int InterfacialState(const int, const int, const int, const int, const int,
        const int [restrict][DIMS], const Node *, const Partition *);
static int GhostState(const int, const int, const int, const int, const int,
        const int [restrict][DIMS], const Node *, const Partition *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * It's necessary to differentiate exterior nodes and interior nodes to avoid
 * incorrect mark of interfacial nodes near global domain boundaries.
 *
 * The identification process should proceed step by step to correctly handle
 * all these relationships and avoid interference between each step,
 * interference may happen if identification processes are crunched together.
 *
 * Moreover, whenever identifying a node, link the node to geometry, which 
 * can be used to access the geometric data for computing required data. 
 * The rational is that don't store every information for each node, but
 * only store necessary information. When need it, access and calculate it.
 */
void ComputeGeometryDomain(Space *space)
{
    ComputeGeometryParameters(space->part.collapsed, &(space->geo));
    InitializeGeometryDomain(space);
    IdentifyGeometryNode(space);
    IdentifyInterfacialNode(space);
    return;
}
static void InitializeGeometryDomain(Space *space)
{
    const Partition *restrict part = &(space->part);
    Node *node = space->node;
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
    Node *node = space->node;
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
 */
static void IdentifyInterfacialNode(Space *space)
{
    const Partition *restrict part = &(space->part);
    Node *node = space->node;
    int idx = 0; /* linear array index math variable */
    const int end = LAYER1 + part->ng * (LAYER2 - LAYER1); /* max search layer is ng + 1 */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                if ((NONE != node[idx].faceID) && (0 == node[idx].geoID)) {
                    /* a newly joined node */
                    FlowReconstruction();
                    node[idx].faceID = NONE; /* reset */
                }
                //if (0 == node[idx].geoID) { /* skip interfacial nodes for main domain */
                //    continue;
                //}
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
        const int path[restrict][DIMS], const Node *node, const Partition *part)
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
        const int path[restrict][DIMS], const Node *node, const Partition *part)
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
    Index idx = 0; /* linear array index math variable */
    Ghost gs = {0.0}; /* data collection for flow reconstruction */
    Real weightSum = 0.0; /* store the sum of weights */
    for (int r = 1; r < space->ng + 2; ++r) {
        for (int k = part->n[PIN][Z][MIN]; k < part->n[PIN][Z][MAX]; ++k) {
            for (int j = part->n[PIN][Y][MIN]; j < part->n[PIN][Y][MAX]; ++j) {
                for (int i = part->n[PIN][X][MIN]; i < part->n[PIN][X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, space);
                    if (GHOST != space->node[idx].type) {
                        continue;
                    }
                    if (r != space->node[idx].layerID) { /* not a ghost node in current layer */
                        continue;
                    }
                    if (model->layers >= r) { /* using the method of image */
                        CalculateGeometryInformation(&gs, k, j, i, geo, space);
                        /* obtain the spatial coordinates of the image point */
                        imageX = info[GSX] + 2 * info[GSDS] * info[GSNX];
                        imageY = info[GSY] + 2 * info[GSDS] * info[GSNY];
                        imageZ = info[GSZ] + 2 * info[GSDS] * info[GSNZ];
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
                        FlowReconstruction(UoImage, imageZ, imageY, imageX, 
                                ComputeK(imageZ, space), ComputeJ(imageY, space), ComputeI(imageX, space),
                                type, UoBC, info, geo, U, space, model, geometry);
                        /* 
                         * Apply the method of image.
                         *  -- reflecting vectors over wall for both slip and noslip, stationary and
                         *     moving conditions is unified by linear interpolation.
                         *  -- scalars are symmetrically reflected between image and ghost.
                         */
                        UoGhost[1] = 2 * UoBC[1] - UoImage[1];
                        UoGhost[2] = 2 * UoBC[2] - UoImage[2];
                        UoGhost[3] = 2 * UoBC[3] - UoImage[3];
                        UoGhost[4] = UoImage[4];
                        UoGhost[5] = UoImage[5];
                    } else { /* using inverse distance weighting instead of method of image */
                        InverseDistanceWeighting(UoGhost, &weightSum, ComputeZ(k, space), ComputeY(j, space), ComputeX(i, space), 
                                k, j, i, 1, type - 1, U, space, model, geometry);
                        /* Normalize the weighted values */
                        NormalizeWeightedValues(UoGhost, weightSum);
                    }
                    UoGhost[0] = UoGhost[4] / (UoGhost[5] * model->gasR); /* compute density */
                    if (!((0 < UoGhost[0]) && (FLT_MAX > UoGhost[0]))) {
                        fprintf(stderr, "rho=%.6g\n", UoGhost[0]);
                        FatalError("illegal density reconstructed, solution diverges...");
                    }
                    ConservativeByPrimitive(U, idx * DIMU, UoGhost, model);
                }
            }
        }
    }
    return;
}
int FlowReconstruction(Real Uo[], const Real z, const Real y, const Real x,
        const int k, const int j, const int i, const int h, 
        Real UoBC[], const Real info[], const Real *geo, const Real *U,
        const Space *space, const Model *model, const Geometry *geometry)
{
    Real weightSum = 0.0; /* store the sum of weights */
    /* pre-estimate step */
    InverseDistanceWeighting(Uo, &weightSum, z, y, x, k, j, i, h, FLUID, U, space, model, geometry);
    /* physical boundary condition enforcement step */
    if (0 < geo[GROUGH]) { /* noslip wall */
        UoBC[1] = geo[GU];
        UoBC[2] = geo[GV];
        UoBC[3] = geo[GW];
    } else { /* slip wall */
        Real nVec[DIMS] = {0.0}; /* normal vector */
        Real taVec[DIMS] = {0.0}; /* tangential vector */
        Real tbVec[DIMS] = {0.0}; /* tangential vector */
        Real rhs[DIMS] = {0.0}; /* right hand side vector */
        OrthogonalSpace(nVec, taVec, tbVec, info);
        rhs[X] = geo[GU] * nVec[X] + geo[GV] * nVec[Y] + geo[GW] * nVec[Z];
        rhs[Y] = (Uo[1] * taVec[X] + Uo[2] * taVec[Y] + Uo[3] * taVec[Z]) / weightSum;
        rhs[Z] = (Uo[1] * tbVec[X] + Uo[2] * tbVec[Y] + Uo[3] * tbVec[Z]) / weightSum;
        UoBC[1] = nVec[X] * rhs[X] + taVec[X] * rhs[Y] + tbVec[X] * rhs[Z];
        UoBC[2] = nVec[Y] * rhs[X] + taVec[Y] * rhs[Y] + tbVec[Y] * rhs[Z];
        UoBC[3] = nVec[Z] * rhs[X] + taVec[Z] * rhs[Y] + tbVec[Z] * rhs[Z];
    }
    UoBC[4] = Uo[4] / weightSum; /* zero gradient pressure */
    if (0 > geo[GT]) { /* adiabatic, dT/dn = 0 */
        UoBC[5] = Uo[5] / weightSum;
    } else { /* otherwise, use specified constant wall temperature, T = Tw */
        UoBC[5] = geo[GT];
    }
    UoBC[0] = UoBC[4] / (UoBC[5] * model->gasR); /* compute density */
    /* correction step by adding the boundary point as a stencil */
    ApplyWeighting(Uo, &weightSum, info[GSDS] * info[GSDS], UoBC, space->tinyL);
    /* Normalize the weighted values */
    NormalizeWeightedValues(Uo, weightSum);
    return 0;
}
static int InverseDistanceWeighting(Real Uo[], Real *weightSum, const Real z, const Real y, const Real x,
        const int k, const int j, const int i, const int h, const int nodeType, const Real *U,
        const Space *space, const Model *model, const Geometry *geometry)
{
    int idxh = 0; /* linear array index math variable */
    int tally = 0; /* stencil count and zero stencil detector */
    int geoID = 0;
    Real Uoh[DIMUo] = {0.0}; /* primitive at target node */
    Real distZ = 0.0;
    Real distY = 0.0;
    Real distX = 0.0;
    Real distance = 0.0;
    for (int dim = 0; dim < DIMUo; ++dim) {
        Uo[dim] = 0.0; /* reset */
    }
    *weightSum = 0.0; /* reset */
    /* 
     * Search around the specified node to find required target nodes as interpolation stencil.
     */
    for (int kh = -h; kh <= h; ++kh) {
        for (int jh = -h; jh <= h; ++jh) {
            for (int ih = -h; ih <= h; ++ih) {
                idxh = IndexMath(k + kh, j + jh, i + ih, space);
                if ((0 > idxh) || (space->nMax <= idxh)) {
                    continue; /* illegal index */
                }
                if (FLUID == nodeType) { /* require fluid nodes */
                    if (FLUID != space->nodeFlag[idxh]) {
                        continue;
                    }
                } else { /* require specified ghost node type */
                    geoID = space->nodeFlag[idxh] - OFFSET - (nodeType - 1) * geometry->totalN;
                    if ((0 > geoID) || (geometry->totalN <= geoID)) { /* not a ghost node with current type */
                        continue;
                    }
                }
                ++tally;
                distZ = ComputeZ(k + kh, space) - z;
                distY = ComputeY(j + jh, space) - y;
                distX = ComputeX(i + ih, space) - x;
                /* use distance square to avoid expensive sqrt */
                distance = distX * distX + distY * distY + distZ * distZ;
                PrimitiveByConservative(Uoh, idxh * DIMU, U, model);
                ApplyWeighting(Uo, weightSum, distance, Uoh, space->tinyL);
            }
        }
    }
    if (0 == tally) {
        fprintf(stderr, "k=%d, j=%d, i=%d\n", k, j, i);
        FatalError("zero stencil encountered...");
    }
    return 0;
}
static int ApplyWeighting(Real Uo[], Real *weightSum, Real distance, const Real Uoh[], const Real tiny)
{
    if (tiny > distance) { /* avoid overflow of too small distance */
        distance = tiny;
    }
    distance = 1 / distance; /* compute weight */
    for (int n = 0; n < DIMUo; ++n) {
        Uo[n] = Uo[n] + Uoh[n] * distance;
    }
    *weightSum = *weightSum + distance; /* accumulate normalizer */
    return 0;
}
/* a good practice: end file with a newline */

