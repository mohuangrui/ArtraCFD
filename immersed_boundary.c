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
    R = 2,
} IBM;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void InitializeGeometryDomain(Space *);
static void IdentifyGeometryNode(Space *);
static void IdentifyInterfacialNode(Space *, const Model *);
static int InterfacialState(const int, const int, const int, const int, const int,
        const int [restrict][DIMS], const Node *const, const Partition *);
static int GhostState(const int, const int, const int, const int, const int,
        const int [restrict][DIMS], const Node *const, const Partition *);
static void ApplyWeighting(const Real [restrict], const Real, Real, 
        Real [restrict], Real [restrict]);
static Real InverseDistanceWeighting(const int, const int [restrict], const Real [restrict],
        const int, const int, const int, const Partition *, const Node *const, const Model *, Real [restrict]);
static void FlowReconstruction(const int, const int [restrict], const Real [restrict], const int,
        const int, const int, const Polyhedron *, const Partition *, const Node *const, const Model *,
        const Real [restrict], const Real [restrict], Real [restrict], Real [restrict]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Mo, H., Lien, F. S., Zhang, F., & Cronin, D. S. (2017). A novel field
 * function for solving complex and dynamic fluid-solid system on Cartesian
 * grid. arXiv:1702.02474.
 *
 * Node mapping should proceed step by step to avoid interference.
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
    Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    int idx = 0; /* linear array index math variable */
    int gid = 0; /* store geometry identifier */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                gid = node[idx].gid;
                if (0 >= gid) {
                    node[idx].gid = 0;
                    node[idx].lid = 0;
                    node[idx].gst = 0;
                    continue;
                }
                poly = geo->poly + gid - 1;
                if (1 == poly->state) {
                    node[idx].lid = 0;
                    node[idx].gst = 0;
                    continue;
                }
                if (0 < node[idx].lid) {
                    node[idx].gid = 0;
                    node[idx].lid = 0;
                    node[idx].gst = 0;
                }
            }
        }
    }
    return;
}
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
    int fid = 0;
    int idx = 0; /* linear array index math variable */
    int box[DIMS][LIMIT] = {{0}}; /* bounding box in node space */
    RealVec p = {0.0}; /* node point */
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        if (1 == poly->state) {
            continue;
        }
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
                    if (0 != node[idx].gid) { /* already classified */
                        continue;
                    }
                    p[X] = PointSpace(i, sMin[X], d[X], ng);
                    p[Y] = PointSpace(j, sMin[Y], d[Y], ng);
                    p[Z] = PointSpace(k, sMin[Z], d[Z], ng);
                    if (0 == poly->faceN) { /* analytical sphere */
                        if (poly->r * poly->r >= Dist2(poly->O, p)) {
                            node[idx].gid = n + 1;
                            node[idx].fid = 0;
                        }
                    } else { /* triangulated polyhedron */
                        if (PointInPolyhedron(p, poly, &fid)) {
                            node[idx].gid = n + 1;
                            node[idx].fid = fid;
                        }
                    }
                }
            }
        }
    }
    return;
}
static void IdentifyInterfacialNode(Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    IntVec n = {0};
    RealVec p = {0.0};
    Real Uo[DIMUo] = {0.0};
    Real weightSum = 0.0;
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                if ((NONE != node[idx].fid) && (0 == node[idx].gid)) {
                    /* a newly joined node */
                    n[X] = i;
                    n[Y] = j;
                    n[Z] = k;
                    p[X] = PointSpace(i, part->domain[X][MIN], part->d[X], part->ng);
                    p[Y] = PointSpace(j, part->domain[Y][MIN], part->d[Y], part->ng);
                    p[Z] = PointSpace(k, part->domain[Z][MIN], part->d[Z], part->ng);
                    weightSum = InverseDistanceWeighting(TO, n, p, R, NONE, 0, part, node, model, Uo);
                    Normalize(DIMUo, weightSum, Uo);
                    Uo[0] = Uo[4] / (Uo[5] * model->gasR); /* compute density */
                    ConservativeByPrimitive(model->gamma, Uo, node[idx].U[TO]);
                    node[idx].fid = NONE; /* reset after correct reconstruction */
                }
                if (0 == node[idx].gid) { /* skip interfacial nodes for main domain */
                    continue;
                }
                /* 
                 * Search neighbours to determine interfacial state. No matter
                 * whether the geometric information of the node is preserved,
                 * the interfacial state of each node should always be reset and
                 * redetermined according to the current domain state.
                 */
                node[idx].lid = InterfacialState(k, j, i, node[idx].gid, part->pathSep[0], part->path, node, part);
                node[idx].gst = 0; /* reset ghost state for potentially uncleaned nodes */
                if ((0 != node[idx].lid) && (0 != node[idx].gid)) { /* an interfacial node may be a ghost node */
                    /* search neighbours to determine ghost state */
                    node[idx].gst = GhostState(k, j, i, 0, part->pathSep[0], part->path, node, part);
                }
            }
        }
    }
    return;
}
static int InterfacialState(const int k, const int j, const int i, const int gid, const int end,
        const int path[restrict][DIMS], const Node *const node, const Partition *part)
{
    /*
     * Search around the specified node to check whether current node
     * is an interfacial node, and return the interfacial state.
     */
    int idx = 0; /* linear array index math variable */
    for (int n = 0; n < end; ++n) {
        idx = IndexNode(k + path[n][Z], j + path[n][Y], i + path[n][X], part->n[Y], part->n[X]);
        if (NONE == node[idx].gid) { /* an exterior node is not valid */
            continue;
        }
        if (gid != node[idx].gid) { /* a heterogeneous node on the path */
            for (int r = 1; r <= part->gl; ++r) {
                if (part->pathSep[r] > n) {
                    return r;
                }
            }
        }
    }
    return 0;
}
static int GhostState(const int k, const int j, const int i, const int gid, const int end,
        const int path[restrict][DIMS], const Node *const node, const Partition *part)
{
    /*
     * Search around the specified node to check whether current node
     * is on numerical boundary, and return the state.
     */
    int idx = 0; /* linear array index math variable */
    for (int n = 0; n < end; ++n) {
        idx = IndexNode(k + path[n][Z], j + path[n][Y], i + path[n][X], part->n[Y], part->n[X]);
        if (NONE == node[idx].gid) { /* an exterior node is not valid */
            continue;
        }
        if (gid == node[idx].gid) { /* a normal computational node on the path */
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
 * Mo, H., Lien, F. S., Zhang, F., & Cronin, D. S. (2016). A novel immersed
 * boundary method for solving flow with arbitrarily irregular and moving
 * geometry. arXiv preprint arXiv:1602.06830.
 */
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
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ValidNodeSpace(NodeSpace(poly->box[s][MIN], sMin[s], dd[s], ng), nMin[s], nMax[s]);
            box[s][MAX] = ValidNodeSpace(NodeSpace(poly->box[s][MAX], sMin[s], dd[s], ng), nMin[s], nMax[s]) + 1;
        }
        for (int r = 1; r <= part->gl; ++r) { /* layer by layer treatment */
            for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
                for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                    for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                        idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                        if ((r != node[idx].gst) || (n + 1 != node[idx].gid)) {
                            continue;
                        }
                        pG[X] = PointSpace(i, sMin[X], d[X], ng);
                        pG[Y] = PointSpace(j, sMin[Y], d[Y], ng);
                        pG[Z] = PointSpace(k, sMin[Z], d[Z], ng);
                        if (model->ibmLayer >= r) { /* immersed boundary treatment */
                            ComputeGeometricData(node[idx].fid, poly, pG, pO, pI, N);
                            nI[X] = NodeSpace(pI[X], sMin[X], dd[X], ng);
                            nI[Y] = NodeSpace(pI[Y], sMin[Y], dd[Y], ng);
                            nI[Z] = NodeSpace(pI[Z], sMin[Z], dd[Z], ng);
                            FlowReconstruction(tn, nI, pI, R, NONE, 0, poly, part, node, model, pO, N, UoO, UoI);
                            MethodOfImage(UoI, UoO, UoG);
                        } else { /* inverse distance weighting */
                            nG[X] = i;
                            nG[Y] = j;
                            nG[Z] = k;
                            weightSum = InverseDistanceWeighting(tn, nG, pG, 1, r - 1, n + 1, part, node, model, UoG);
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
    UoG[1] = UoO[1] + UoO[1] - UoI[1];
    UoG[2] = UoO[2] + UoO[2] - UoI[2];
    UoG[3] = UoO[3] + UoO[3] - UoI[3];
    UoG[4] = UoI[4];
    UoG[5] = UoI[5];
    return;
}
void ComputeGeometricData(const int fid, const Polyhedron *poly, const Real pG[restrict],
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
        ComputeIntersection(pG, fid, poly, pO, N);
    }
    pI[X] = pO[X] + pO[X] - pG[X];
    pI[Y] = pO[Y] + pO[Y] - pG[Y];
    pI[Z] = pO[Z] + pO[Z] - pG[Z];
    return;
}
static void FlowReconstruction(const int tn, const int n[restrict], const Real p[restrict], const int h,
        const int type, const int gid, const Polyhedron *poly, const Partition *part, const Node *const node, 
        const Model *model, const Real pO[restrict], const Real N[restrict], Real UoO[restrict], Real UoI[restrict])
{
    const Real zero = 0.0;
    const Real one = 1.0;
    /* pre-estimate step */
    Real weightSum = InverseDistanceWeighting(tn, n, p, h, type, gid, part, node, model, UoI);
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
        const RealVec VI = {UoI[1] * weight, UoI[2] * weight, UoI[3] * weight};
        RealVec Ta = {zero}; /* tangential vector */
        RealVec Tb = {zero}; /* tangential vector */
        Real RHS[DIMS] = {zero}; /* right hand side vector */
        OrthogonalSpace(N, Ta, Tb);
        RHS[X] = Dot(Vs, N);
        RHS[Y] = Dot(VI, Ta);
        RHS[Z] = Dot(VI, Tb);
        UoO[1] = N[X] * RHS[X] + Ta[X] * RHS[Y] + Tb[X] * RHS[Z];
        UoO[2] = N[Y] * RHS[X] + Ta[Y] * RHS[Y] + Tb[Y] * RHS[Z];
        UoO[3] = N[Z] * RHS[X] + Ta[Z] * RHS[Y] + Tb[Z] * RHS[Z];
    }
    UoO[4] = UoI[4] * weight;
    if (zero > poly->T) { /* adiabatic, dT/dn = 0 */
        UoO[5] = UoI[5] * weight;
    } else { /* otherwise, use specified constant wall temperature, T = Tw */
        UoO[5] = poly->T;
    }
    /* correction step by adding the boundary point as a stencil */
    ApplyWeighting(UoO, part->tinyL, Dist2(p, pO), &weightSum, UoI);
    /* Normalize the weighted values */
    Normalize(DIMUo, weightSum, UoI);
    return;
}
static Real InverseDistanceWeighting(const int tn, const int n[restrict], const Real p[restrict], 
        const int h, const int type, const int gid, const Partition *part, 
        const Node *const node, const Model *model, Real Uo[restrict])
{
    int idx = 0; /* linear array index math variable */
    const int idxMax = part->n[X] * part->n[Y] * part->n[Z];
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const int ng = part->ng;
    Real Uoh[DIMUo] = {0.0}; /* primitive at neighbouring node */
    RealVec ph = {0.0}; /* neighbouring point */
    Real weightSum = 0.0;
    memset(Uo, 0, DIMUo * sizeof(*Uo));
    for (int r = h, tally = 0; 0 == tally; ++r) {
        for (int kh = -r; kh <= r; ++kh) {
            for (int jh = -r; jh <= r; ++jh) {
                for (int ih = -r; ih <= r; ++ih) {
                    idx = IndexNode(n[Z] + kh, n[Y] + jh, n[X] + ih, part->n[Y], part->n[X]);
                    if ((0 > idx) || (idxMax <= idx)) { /* illegal index */
                        continue;
                    }
                    if (gid != node[idx].gid) {
                        continue;
                    }
                    if (0 == gid) { /* require normal node type */
                        if (type != node[idx].fid) {
                            continue;
                        }
                    } else { /* require specified ghost node type */
                        if (type != node[idx].gst) { /* not a ghost node with current type */
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

