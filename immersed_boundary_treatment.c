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
#include "immersed_boundary_treatment.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <float.h> /* size of floating point values */
#include "computational_geometry.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int InitializeGeometryDomain(Space *, const Partition *);
static int IdentifySolidNodes(Space *, const Partition *, const Geometry *);
static int IdentifyGhostNodes(Space *, const Partition *, const Geometry *);
static int SearchFluidNodes(const int, const int, const int, const int, const Space *);
static int InverseDistanceWeighting(Real [], Real *, const Real, const Real, const Real,
        const int, const int, const int, const int, const int, const Real *,
        const Space *, const Model *, const Geometry *);
static int ApplyWeighting(Real [], Real *, Real, const Real [], const Real);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function identify the type of each node: 
 *
 * Procedures are:
 * -- initialize node flag of global boundary and exterior nodes to exterior
 *    type, and inner nodes to fluid type;
 * -- identify all inner nodes that are in geometry as solid node;
 * -- identify ghost nodes by the node type of neighbours of solid nodes;
 *
 * Ghost nodes are solid nodes that locate on numerical boundaries.
 * It's necessary to differ global boundary nodes and inner nodes to avoid 
 * incorrect mark of ghost nodes nearby global domain boundaries.
 *
 * The identification process should proceed step by step to correctly handle
 * all these relationships and avoid interference between each step,
 * interference may happen if identification processes are crunched together.
 *
 * Moreover, whenever identifying a solid node or ghost node, store its
 * corresponding geometry information by linking the node flag to the
 * geometry ID, which will be accessed to calculate other informations. 
 * The rational is that don't store every information for each ghost node, but
 * only store necessary information. When need it, access and calculate it.
 */
int ComputeGeometryDomain(Space *space, const Partition *part, const Geometry *geo)
{
    InitializeGeometryDomain(space, part);
    IdentifySolidNodes(space, part, geo);
    IdentifyGhostNodes(space, part, geo);
    return 0;
}
static int InitializeGeometryDomain(Space *space, const Partition *part)
{
    Index idx = 0; /* linear array index math variable */
    /* initialize inner nodes to fluid type, others to exterior nodes type. */
    for (int k = 0; k < space->n[Z]; ++k) {
        for (int j = 0; j < space->n[Y]; ++j) {
            for (int i = 0; i < space->n[X]; ++i) {
                idx = IndexNode(k, j, i, space);
                if ((part->n[PIN][Z][MIN] <= k) && (part->n[PIN][Z][MAX] > k) &&
                        (part->n[PIN][Y][MIN] <= j) && (part->n[PIN][Y][MAX] > j) &&
                        (part->n[PIN][X][MIN] <= i) && (part->n[PIN][X][MAX] > i)) {
                    space->node[idx].type = FLUID;
                } else {
                    space->node[idx].type = EXTERIOR;
                }
            }
        }
    }
    return 0;
}
/*
 * When locate solid nodes, there are two approaches available. One is search
 * over each node and verify each node regarding to all the geometries; another
 * is search each geometry and find all the nodes inside current geometry.
 * The second method is adopted here for performance reason, although it's much
 * more complicated than the first one.
 * It is efficient to only test points that are inside the bounding box or
 * sphere of a large polyhedron. Be cautious with the validity of any calculated
 * index. It's extremely necessary to adjust the index into the valid region or 
 * check the validity of the index to avoid index exceed array bound limits.
 */
static int IdentifySolidNodes(Space *space, const Partition *part, const Geometry *geo)
{
    Index idx = 0; /* linear array index math variable */
    int box[DIMS][LIMIT] = {{0}}; /* range box in node space */
    for (int m = 0; m < geo->totalM; ++m) {
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ValidNodeSpace(NodeSpace(geo->list[m].box[s][MIN], s, space), s, part);
            box[s][MAX] = ValidNodeSpace(NodeSpace(geo->list[m].box[s][MAX], s, space), s, part) + 1;
        }
        /* find nodes in geometry, then flag as solid and link to geometry. */
        for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
            for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                    if (0 > PointInPolyhedron(k, j, i, geo->list + m, space)) {
                        idx = IndexNode(k, j, i, space);
                        space->node[idx].type = SOLID;
                        space->node[idx].geoID = m;
                    }
                }
            }
        }
    }
    return 0;
}
static int IdentifyGhostNodes(Space *space, const Partition *part, const Geometry *geo)
{
    Index idx = 0; /* linear array index math variable */
    /* find solid nodes at numerical boundary, then flag as ghost while preserving link. */
    for (int k = part->n[PIN][Z][MIN]; k < part->n[PIN][Z][MAX]; ++k) {
        for (int j = part->n[PIN][Y][MIN]; j < part->n[PIN][Y][MAX]; ++j) {
            for (int i = part->n[PIN][X][MIN]; i < part->n[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, space);
                if (SOLID != space->node[idx].type) {
                    continue;
                }
                /* search neighbours of node(k, j, i) to check fluid node */
                for (int r = 1; r < space->ng + 2; ++r) { /* max search range should be ng + 1 */
                    /* if rth neighbours has fluid, then it's rth type ghost node */
                    if (0 == SearchFluidNodes(k, j, i, r, space)) {
                        space->node[idx].type = GHOST;
                        space->node[idx].layerID = r; /* save ghost layer identifier */
                        break; /* exit search loop once identified successfully */
                    }
                }
            }
        }
    }
    return 0;
}
static int SearchFluidNodes(const int k, const int j, const int i, 
        const int h, const Space *space)
{
    /*
     * Search around the specified node and return the information of whether 
     * current node has fluid neighbours at the specified coordinate range.
     */
    const Index idxW = IndexNode(k, j, i - h, space);
    const Index idxE = IndexNode(k, j, i + h, space);
    const Index idxS = IndexNode(k, j - h, i, space);
    const Index idxN = IndexNode(k, j + h, i, space);
    const Index idxF = IndexNode(k - h, j, i, space);
    const Index idxB = IndexNode(k + h, j, i, space);
    return (space->node[idxW].type * space->node[idxE].type * 
            space->node[idxS].type * space->node[idxN].type * 
            space->node[idxF].type * space->node[idxB].type);
}
/*
 * Numerical boundary treatments for ghost nodes.
 */
int BoundaryTreatmentsGCIBM(Real *U, const Space *space, const Model *model,
        const Partition *part, const Geometry *geometry)
{
    Real UoGhost[DIMUo] = {0.0}; /* reconstructed primitives at ghost node */
    Real UoImage[DIMUo] = {0.0}; /* reconstructed primitives at image point */
    Real UoBC[DIMUo] = {0.0}; /* physical primitives at boundary point */
    Real info[INFOGHOST] = {0.0}; /* store calculated geometry information */
    Real weightSum = 0.0; /* store the sum of weights */
    Real imageZ = 0.0;
    Real imageY = 0.0;
    Real imageX = 0.0;
    int idx = 0; /* linear array index math variable */
    int geoID = 0; /* geometry id */
    Real *geo = NULL;
    for (int type = 1; type < space->ng + 2; ++type) {
        for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
            for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
                for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                    idx = IndexMath(k, j, i, space);
                    if (OFFSET > space->nodeFlag[idx]) { /* it's not a ghost */
                        continue;
                    }
                    geoID = space->nodeFlag[idx] - OFFSET - (type - 1) * geometry->totalN;
                    if ((0 > geoID) || (geometry->totalN <= geoID)) { /* not a ghost node with current type */
                        continue;
                    }
                    if (model->layers >= type) { /* using method of image */
                        geo = IndexGeometry(geoID, geometry);
                        CalculateGeometryInformation(info, k, j, i, geo, space);
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
    return 0;
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
Real InGeometry(const int k, const int j, const int i, const Real *geo, const Space *space)
{
    /* x, y, z distance to center */
    const Real lX = ComputeX(i, space) - geo[GX];
    const Real lY = ComputeY(j, space) - geo[GY];
    const Real lZ = ComputeZ(k, space) - geo[GZ];
    return (lX * lX + lY * lY + lZ * lZ - geo[GR] * geo[GR]);
}
int CalculateGeometryInformation(Real info[], const int k, const int j, const int i, 
        const Real *geo, const Space *space)
{
    /* x, y, z coordinates */
    info[GSX] = ComputeX(i, space);
    info[GSY] = ComputeY(j, space);
    info[GSZ] = ComputeZ(k, space);
    /* temporary store the x, y, z distance to geometry center */
    info[GSNX] = info[GSX] - geo[GX];
    info[GSNY] = info[GSY] - geo[GY];
    info[GSNZ] = info[GSZ] - geo[GZ];
    /* distance to center */
    info[GSDC] = sqrt(info[GSNX] * info[GSNX] + info[GSNY] * info[GSNY] + info[GSNZ] * info[GSNZ]);
    /* distance to surface */
    info[GSDS] = geo[GR] - info[GSDC];
    /* x, y, z normal vector components */
    info[GSNX] = info[GSNX] / info[GSDC];
    info[GSNY] = info[GSNY] / info[GSDC];
    info[GSNZ] = info[GSNZ] / info[GSDC];
    return 0;
}
/* a good practice: end file with a newline */

