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
#include <stdlib.h> /* mathematical functions on integers */
#include <float.h> /* size of floating point values */
#include <string.h> /* manipulating strings */
#include "computational_geometry.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    R = 2, /* domain of dependence radius */
    INTERL = 0, /* interfacial layer state */
    INTERG = 1, /* ghost layer state */
    TYPED = -1, /* domain as key reconstruction state */
    TYPEF = -2, /* face as key reconstruction state */
    TYPEL = -3, /* layer as key reconstruction state */
} IbmConst;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void InitializeGeometricField(Space *);
static void SetDomainField(Space *);
static void SetInterfacialField(Space *, const Model *);
static int GetInterState(const int, const int, const int, const int, const int,
        const int, const int [restrict][DIMS], const Node *const, const Partition *const);
static void ApplyWeighting(const Real [restrict], const Real, Real,
        Real [restrict], Real [restrict]);
static Real InverseDistanceWeighting(const int, const int [restrict],
        const Real [restrict], const int, const int, const int, const Partition *const,
        const Node *const, const Model *, Real [restrict]);
static void ReconstructFlow(const int, const int [restrict], const Real [restrict],
        const int, const int, const int, const Polyhedron *, const Partition *const,
        const Node *const, const Model *, const Real [restrict], const Real [restrict],
        Real [restrict], Real [restrict]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Mo, H., Lien, F. S., Zhang, F., & Cronin, D. S. (2017). A novel field
 * function for solving complex and dynamic fluid-solid system on Cartesian
 * grid. arXiv:1702.02474.
 *
 * Node mapping should proceed step by step to avoid interference.
 * Once a node is mapped, it is linked to geometry. This link is used to
 * access the geometry for computing future required geometric data.
 * The rational is that do not store every information for each node, but
 * only store necessary links. When need it, access and calculate it.
 */
void ComputeGeometricField(Space *space, const Model *model)
{
    InitializeGeometricField(space);
    SetDomainField(space);
    SetInterfacialField(space, model);
    return;
}
static void InitializeGeometricField(Space *space)
{
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    const Geometry *const geo = &(space->geo);
    const Polyhedron *poly = NULL;
    int idx = 0; /* linear array index math variable */
    int gid = 0; /* store geometry identifier */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                gid = node[idx].did;
                node[idx].gst = node[idx].did; /* preserve domain field */
                if (0 == gid) {
                    node[idx].fid = 0; /* remove passe domain change mark */
                    continue; /* skip non-polyhedron nodes */
                }
                /* the rest is to treat polyhedron nodes */
                poly = geo->poly + gid - 1;
                if (1 == poly->state) {
                    continue; /* keep domain field for nodes in stationary polyhedron */
                }
                /*
                 * The rest is to treat nodes in non-stationary polyhedrons. Due
                 * to the restricted motion, can only reset interfacial nodes for
                 * remapping while keeping non-interfacial nodes to reduce cost.
                 * When polyhedrons move, the previous nth layer may become a
                 * (n-1)th layer, therefore, need to reset gl+1 layers to
                 * ensure the closest face id information of all the future
                 * gl interfacial nodes are updated. However, if only need to
                 * update the closest face id information for the future gl-1
                 * layers, can only reset gl interfacial layers.
                 */
                if (0 < node[idx].lid) {
                    node[idx].did = 0;
                }
            }
        }
    }
    return;
}
/*
 * In domain-node mapping, there are two approaches available. One is loop
 * over each node to verify each node regarding all the geometries; another
 * is loop over each geometry to find all the nodes inside current geometry.
 * The second method is adopted here for performance reason.
 *
 * To best utilize the convergence property of immersed boundary treatment,
 * points either in or on geometry should be classified into the geometry.
 *
 * For a large data set, an extra data preprocessing with spatial subdivision
 * is required to ensure minimal test in addition to the bounding container
 * method. Spatial subdivision is to provide internal resolution for the
 * polyhedron for fast inclusion determination.
 */
static void SetDomainField(Space *space)
{
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    const Geometry *const geo = &(space->geo);
    const IntVec nMin = {part->ns[PIN][X][MIN], part->ns[PIN][Y][MIN], part->ns[PIN][Z][MIN]};
    const IntVec nMax = {part->ns[PIN][X][MAX], part->ns[PIN][Y][MAX], part->ns[PIN][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const IntVec ng = {part->ng[X], part->ng[Y], part->ng[Z]};
    const Polyhedron *poly = NULL;
    int box[DIMS][LIMIT] = {{0}}; /* bounding box in node space */
    int fid = 0; /* store face link */
    int idx = 0; /* linear array index math variable */
    RealVec p = {0.0}; /* node point */
    /* overlapping geometries introduce loop-carried dependence for node mapping */
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        if (1 == poly->state) {
            continue;
        }
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ConfineSpace(MapNode(poly->box[s][MIN], sMin[s], dd[s], ng[s]), nMin[s], nMax[s]);
            box[s][MAX] = ConfineSpace(MapNode(poly->box[s][MAX], sMin[s], dd[s], ng[s]), nMin[s], nMax[s]) + 1;
        }
        /* find nodes in geometry, then flag and link to geometry */
        for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
            for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    if (0 != node[idx].did) { /* already classified */
                        continue;
                    }
                    p[X] = MapPoint(i, sMin[X], d[X], ng[X]);
                    p[Y] = MapPoint(j, sMin[Y], d[Y], ng[Y]);
                    p[Z] = MapPoint(k, sMin[Z], d[Z], ng[Z]);
                    if (0 >= poly->faceN) { /* analytical polyhedron */
                        if (poly->r * poly->r >= Dist2(poly->O, p)) {
                            node[idx].did = n + 1;
                            node[idx].fid = 0;
                        }
                    } else { /* triangulated polyhedron */
                        if (PointInPolyhedron(p, poly, &fid)) {
                            node[idx].did = n + 1;
                            node[idx].fid = fid;
                        }
                    }
                }
            }
        }
    }
    return;
}
static void SetInterfacialField(Space *space, const Model *model)
{
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    const int sd = 0; /* solution domain */
    IntVec n = {0}; /* current node */
    RealVec p = {0.0}; /* node point */
    Real Uo[DIMUo] = {0.0};
    Real weightSum = 0.0;
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                /* reconstruct newly joined node for the solution domain */
                if ((node[idx].gst != node[idx].did) && (sd == node[idx].did)) {
                    /* a newly joined solution domain node */
                    n[X] = i; n[Y] = j; n[Z] = k;
                    p[X] = MapPoint(i, part->domain[X][MIN], part->d[X], part->ng[X]);
                    p[Y] = MapPoint(j, part->domain[Y][MIN], part->d[Y], part->ng[Y]);
                    p[Z] = MapPoint(k, part->domain[Z][MIN], part->d[Z], part->ng[Z]);
                    weightSum = InverseDistanceWeighting(TO, n, p, R, TYPEF, node[idx].did, part, node, model, Uo);
                    Normalize(DIMUo, weightSum, Uo);
                    Uo[0] = Uo[4] / (Uo[5] * model->gasR); /* compute density */
                    MapConservative(model->gamma, Uo, node[idx].U[TO]);
                    node[idx].fid = NONE; /* set domain change mark to avoid reconstruction interference */
                }
                /* reset interfacial state */
                node[idx].lid = 0;
                node[idx].gst = 0;
                /* search neighbours to determine the current interfacial state */
                if (sd == node[idx].did) { /* skip interfacial nodes for main domain */
                    continue;
                }
                node[idx].lid = GetInterState(INTERL, k, j, i, node[idx].did, part->pathSep[0], part->path, node, part);
                if ((0 < node[idx].lid) && (sd != node[idx].did)) { /* ghost node is a subset of interfacial node */
                    node[idx].gst = GetInterState(INTERG, k, j, i, sd, part->pathSep[0], part->path, node, part);
                }
            }
        }
    }
    return;
}
static int GetInterState(const int sid, const int k, const int j, const int i, const int did,
        const int end, const int path[restrict][DIMS], const Node *const node, const Partition *const part)
{
    /* search around the specified node to check interfacial state */
    int idx = 0; /* linear array index math variable */
    int ih = 0, jh = 0, kh = 0; /* neighbouring node */
    int flag = 0; /* control flag */
    for (int n = 0; n < end; ++n) {
        kh = k + path[n][Z];
        jh = j + path[n][Y];
        ih = i + path[n][X];
        if (!InPartBox(kh, jh, ih, part->ns[PIN])) {
            continue;
        }
        idx = IndexNode(kh, jh, ih, part->n[Y], part->n[X]);
        switch (sid) {
            case INTERL:
                if (did != node[idx].did) { /* a heterogeneous node on the path */
                    flag = 1;
                }
                break;
            case INTERG:
                if (did == node[idx].did) { /* a computational node on the path */
                    flag = 1;
                }
                break;
            default:
                break;
        }
        if (1 == flag) { /* return the layer number */
            for (int r = 1; r <= part->gl; ++r) {
                if (part->pathSep[r] > n) {
                    return r;
                }
            }
        }
    }
    return 0;
}
/*
 * Mo, H., Lien, F.S., Zhang, F. and Cronin, D.S., 2016. A sharp interface
 * immersed boundary method for solving flow with arbitrarily irregular and
 * changing geometry. arXiv:1602.06830.
 */
void TreatImmersedBoundary(const int tn, Space *space, const Model *model)
{
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    const Geometry *const geo = &(space->geo);
    const IntVec nMin = {part->ns[PIN][X][MIN], part->ns[PIN][Y][MIN], part->ns[PIN][Z][MIN]};
    const IntVec nMax = {part->ns[PIN][X][MAX], part->ns[PIN][Y][MAX], part->ns[PIN][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const IntVec ng = {part->ng[X], part->ng[Y], part->ng[Z]};
    const Polyhedron *poly = NULL;
    int idx = 0; /* linear array index math variable */
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
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ConfineSpace(MapNode(poly->box[s][MIN], sMin[s], dd[s], ng[s]), nMin[s], nMax[s]);
            box[s][MAX] = ConfineSpace(MapNode(poly->box[s][MAX], sMin[s], dd[s], ng[s]), nMin[s], nMax[s]) + 1;
        }
        /* treat ghost nodes */
        for (int r = 1; r <= part->gl; ++r) { /* layer by layer treatment */
            for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
                for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                    for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                        idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                        if ((r != node[idx].gst) || (n + 1 != node[idx].did)) {
                            continue;
                        }
                        pG[X] = MapPoint(i, sMin[X], d[X], ng[X]);
                        pG[Y] = MapPoint(j, sMin[Y], d[Y], ng[Y]);
                        pG[Z] = MapPoint(k, sMin[Z], d[Z], ng[Z]);
                        if (model->ibmLayer >= r) { /* immersed boundary treatment */
                            ComputeGeometricData(pG, node[idx].fid, poly, pO, pI, N);
                            nI[X] = MapNode(pI[X], sMin[X], dd[X], ng[X]);
                            nI[Y] = MapNode(pI[Y], sMin[Y], dd[Y], ng[Y]);
                            nI[Z] = MapNode(pI[Z], sMin[Z], dd[Z], ng[Z]);
                            /*
                             * When extremely strong discontinuities exist in the
                             * domain of dependence of inverse distance weighting,
                             * WENO's idea may be adopted to avoid discontinuous
                             * stencils and to only use smooth stencils. However,
                             * the algorithm will be too complex.
                             */
                            ReconstructFlow(tn, nI, pI, R, TYPED, 0, poly, part, node, model, pO, N, UoO, UoI);
                            DoMethodOfImage(UoI, UoO, UoG);
                        } else { /* inverse distance weighting */
                            nG[X] = i; nG[Y] = j; nG[Z] = k;
                            weightSum = InverseDistanceWeighting(tn, nG, pG, 1, r - 1, n + 1, part, node, model, UoG);
                            Normalize(DIMUo, weightSum, UoG);
                        }
                        UoG[0] = UoG[4] / (UoG[5] * model->gasR); /* compute density */
                        MapConservative(model->gamma, UoG, node[idx].U[tn]);
                    }
                }
            }
        }
    }
    return;
}
void DoMethodOfImage(const Real UoI[restrict], const Real UoO[restrict], Real UoG[restrict])
{
    /*
     * Apply the method of image.
     *  -- reflecting vectors over wall for slip and noslip as well as stationary
     *     and moving conditions is unified by linear interpolation.
     *  -- scalars are symmetrically reflected between image and ghost.
     *  -- other scalars are determined by equation of state.
     */
    UoG[1] = UoO[1] + UoO[1] - UoI[1];
    UoG[2] = UoO[2] + UoO[2] - UoI[2];
    UoG[3] = UoO[3] + UoO[3] - UoI[3];
    UoG[4] = UoI[4];
    UoG[5] = UoI[5];
    return;
}
static void ReconstructFlow(const int tn, const int n[restrict], const Real p[restrict],
        const int h, const int type, const int did, const Polyhedron *poly, const Partition *const part,
        const Node *const node, const Model *model, const Real pO[restrict], const Real N[restrict],
        Real UoO[restrict], Real Uo[restrict])
{
    const Real zero = 0.0;
    const Real one = 1.0;
    /* pre-estimate step */
    Real weightSum = InverseDistanceWeighting(tn, n, p, h, type, did, part, node, model, Uo);
    const Real weight = one / weightSum;
    /* physical boundary condition enforcement step */
    RealVec Vs = {zero}; /* general motion of boundary point */
    /* Vs = Vcentroid + W x r */
    const RealVec r = {pO[X] - poly->O[X], pO[Y] - poly->O[Y], pO[Z] - poly->O[Z]};
    Cross(poly->W[TO], r, Vs); /* relative motion in translating coordinate system */
    Vs[X] = poly->V[TO][X] + Vs[X];
    Vs[Y] = poly->V[TO][Y] + Vs[Y];
    Vs[Z] = poly->V[TO][Z] + Vs[Z];
    if (zero < poly->cf) { /* noslip wall */
        UoO[1] = Vs[X];
        UoO[2] = Vs[Y];
        UoO[3] = Vs[Z];
    } else { /* slip wall */
        const RealVec V = {Uo[1] * weight, Uo[2] * weight, Uo[3] * weight};
        RealVec Ta = {zero}; /* tangential vector */
        RealVec Tb = {zero}; /* tangential vector */
        Real RHS[DIMS] = {zero}; /* right hand side vector */
        OrthogonalSpace(N, Ta, Tb);
        RHS[X] = Dot(Vs, N);
        RHS[Y] = Dot(V, Ta);
        RHS[Z] = Dot(V, Tb);
        UoO[1] = N[X] * RHS[X] + Ta[X] * RHS[Y] + Tb[X] * RHS[Z];
        UoO[2] = N[Y] * RHS[X] + Ta[Y] * RHS[Y] + Tb[Y] * RHS[Z];
        UoO[3] = N[Z] * RHS[X] + Ta[Z] * RHS[Y] + Tb[Z] * RHS[Z];
    }
    /*
     * dp/dn = rho_f * v_t^2 / R - rho_f * a_s
     * v_t = (v_o - v_s) - (v_o - vs) * n: relative tangential velocity
     * v_o velocity of boundary point; v_s: velocity of solid surface
     * R: radius of curvature;
     * a_s = at + ar x r + w x (w x r): acceleration of solid surface
     * currently, the effect of this condition is found very limited.
     * Therefore, boundary layer pressure condition dp/dn = 0 is used.
     */
    UoO[4] = Uo[4] * weight;
    if (zero > poly->T) { /* adiabatic, dT/dn = 0 */
        UoO[5] = Uo[5] * weight;
    } else { /* otherwise, use specified constant wall temperature, T = Tw */
        UoO[5] = poly->T;
    }
    /* correction step by adding the boundary point as a stencil */
    ApplyWeighting(UoO, part->tinyL, Dist2(p, pO), &weightSum, Uo);
    /* Normalize the weighted values */
    Normalize(DIMUo, weightSum, Uo);
    return;
}
static Real InverseDistanceWeighting(const int tn, const int n[restrict], const Real p[restrict],
        const int h, const int type, const int did, const Partition *const part,
        const Node *const node, const Model *model, Real Uo[restrict])
{
    int idx = 0; /* linear array index math variable */
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const IntVec ng = {part->ng[X], part->ng[Y], part->ng[Z]};
    Real Uoh[DIMUo] = {0.0}; /* primitive at neighbouring node */
    RealVec ph = {0.0}; /* neighbouring point */
    IntVec nh = {0}; /* neighbouring node */
    Real weightSum = 0.0;
    memset(Uo, 0, DIMUo * sizeof(*Uo));
    /*
     * Search nodes with required "type" in the domain specified by the center node
     * "n" and initial range "h" as interpolation stencils for the interpolated point "p".
     * To preserve symmetry, the search range in each direction should be symmetric, and
     * the search operator on each direction index should be symmetric. In addition,
     * any temporal priority should be strictly prevented in treating the solution nodes.
     */
    for (int r = h, tally = 0; 0 == tally; ++r) {
        for (int kh = -r; kh <= r; ++kh) {
            for (int jh = -r; jh <= r; ++jh) {
                for (int ih = -r; ih <= r; ++ih) {
                    nh[X] = n[X] + ih;
                    nh[Y] = n[Y] + jh;
                    nh[Z] = n[Z] + kh;
                    if (!InPartBox(nh[Z], nh[Y], nh[X], part->ns[PIN])) {
                        continue;
                    }
                    idx = IndexNode(nh[Z], nh[Y], nh[X], part->n[Y], part->n[X]);
                    /* be aware of the validity of ih = jh = kh = 0 */
                    if (did != node[idx].did) {
                        continue; /* skip node not in target domain */
                    }
                    switch (type) {
                        case TYPED: /* use node in target domain */
                            break;
                        case TYPEF: /* use original node in target domain to avoid priority */
                            if ((did != node[idx].gst) || (0 > node[idx].fid)) {
                                continue; /* skip changed node either reconstructed or not */
                            }
                            break;
                        default: /* use node in target domain layer */
                            if (type != node[idx].gst) {
                                continue;
                            }
                            break;
                    }
                    ++tally;
                    ph[X] = MapPoint(nh[X], sMin[X], d[X], ng[X]);
                    ph[Y] = MapPoint(nh[Y], sMin[Y], d[Y], ng[Y]);
                    ph[Z] = MapPoint(nh[Z], sMin[Z], d[Z], ng[Z]);
                    MapPrimitive(model->gamma, model->gasR, node[idx].U[tn], Uoh);
                    ApplyWeighting(Uoh, part->tinyL, Dist2(p, ph), &weightSum, Uo);
                }
            }
        }
    }
    return weightSum;
}
static void ApplyWeighting(const Real Uoh[restrict], const Real tiny, Real weight,
        Real weightSum[restrict], Real Uo[restrict])
{
    const Real one = 1.0;
    if (tiny > weight) { /* avoid overflow of too small weight */
        weight = tiny;
    }
    weight = one / weight; /* compute weight */
    for (int n = 0; n < DIMUo; ++n) {
        Uo[n] = Uo[n] + Uoh[n] * weight;
    }
    *weightSum = *weightSum + weight; /* accumulate normalizer */
    return;
}
/* a good practice: end file with a newline */

