/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) 2014-2018 Huangrui Mo <huangrui.mo@gmail.com>              *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "gcibm.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <float.h> /* size of floating point values */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int InitializeDomainGeometry(Space *, const Partition *);
static int IdentifySolidNodes(Space *, const Partition *, const Geometry *);
static int IdentifyGhostNodes(Space *, const Partition *, const Geometry *);
static int SearchFluidNodes(const int, const int, const int, const int, const Space *);
static int ApplyWeighting(Real [], Real *, Real, const Real [], const Real);
static int OrthogonalSpace(Real *, Real *, Real *, const Real *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * These functions identify the type of each node: 
 *
 * Procedures are:
 * -- initialize node flag of boundary and exterior nodes to exterior type,
 *    and inner nodes to fluid type;
 * -- identify all inner nodes that are in solid geometry as solid node;
 * -- identify ghost nodes according to the node type of its neighbours;
 *
 * Ghost nodes are solid nodes that locate on numerical boundaries.
 * It's necessary to differ boundary nodes and inner nodes to avoid 
 * incorrect mark of ghost nodes nearby domain boundaries.
 * It's necessary to differ different types of ghost node according to
 * which neighbour of them is a fluid node. The fist type of ghost node 
 * is special because the computation of boundary forces are only based on
 * the first nearest solid layer of boundary.
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
int ComputeDomainGeometryGCIBM(Space *space, const Partition *part, const Geometry *geometry)
{
    InitializeDomainGeometry(space, part);
    IdentifySolidNodes(space, part, geometry);
    IdentifyGhostNodes(space, part, geometry);
    return 0;
}
static int InitializeDomainGeometry(Space *space, const Partition *part)
{
    int idx = 0; /* linear array index math variable */
    /* initialize inner nodes to fluid type, others to exterior nodes type. */
    for (int k = 0; k < space->kMax; ++k) {
        for (int j = 0; j < space->jMax; ++j) {
            for (int i = 0; i < space->iMax; ++i) {
                idx = IndexMath(k, j, i, space);
                if ((part->kSub[0] <= k) && (part->kSup[0] > k) &&
                        (part->jSub[0] <= j) && (part->jSup[0] > j) &&
                        (part->iSub[0] <= i) && (part->iSup[0] > i)) {
                    space->nodeFlag[idx] = FLUID;
                } else {
                    space->nodeFlag[idx] = EXTERIOR;
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
 * Be cautious with the validity of any calculated index. It's extremely
 * necessary to adjust the index into the valid region or check the
 * validity of the index to avoid index exceed array bound limits and 
 * mysterious bugs.
 */
static int IdentifySolidNodes(Space *space, const Partition *part, const Geometry *geometry)
{
    Real *geo = NULL;
    Real safetyCoe = 1.2; /* zoom the search range */
    int idx = 0; /* linear array index math variable */
    int centerK = 0;
    int centerJ = 0;
    int centerI = 0;
    int rangeK = 0;
    int rangeJ = 0;
    int rangeI = 0;
    int kSub = 0;
    int kSup = 0;
    int jSub = 0;
    int jSup = 0;
    int iSub = 0;
    int iSup = 0;
    for (int geoCount = 0; geoCount < geometry->totalN; ++geoCount) {
        geo = IndexGeometry(geoCount, geometry);
        centerK = ComputeK(geo[GZ], space);
        centerJ = ComputeJ(geo[GY], space);
        centerI = ComputeI(geo[GX], space);
        rangeK = (int)(safetyCoe * geo[GR] * space->ddz);
        rangeJ = (int)(safetyCoe * geo[GR] * space->ddy);
        rangeI = (int)(safetyCoe * geo[GR] * space->ddx);
        /* determine search range according to valid region */
        kSub = ValidRegionK(centerK - rangeK, part);
        kSup = ValidRegionK(centerK + rangeK, part) + 1;
        jSub = ValidRegionJ(centerJ - rangeJ, part);
        jSup = ValidRegionJ(centerJ + rangeJ, part) + 1;
        iSub = ValidRegionI(centerI - rangeI, part);
        iSup = ValidRegionI(centerI + rangeI, part) + 1;
        for (int k = kSub; k < kSup; ++k) {
            for (int j = jSub; j < jSup; ++j) {
                for (int i = iSub; i < iSup; ++i) {
                    if (0 > InGeometry(k, j, i, geo, space)) {
                        idx = IndexMath(k, j, i, space);
                        space->nodeFlag[idx] = -OFFSET - geoCount; /* geometry are linked */
                    }
                }
            }
        }
    }
    return 0;
}
static int IdentifyGhostNodes(Space *space, const Partition *part, const Geometry *geometry)
{
    int idx = 0; /* linear array index math variable */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (-OFFSET < space->nodeFlag[idx]) { /* it's not solid node */
                    continue;
                }
                /* if nth neighbour is fluid, then it's nth type ghost node */
                for (int type = 1; type < space->ng + 2; ++type) { /* max search range should be ng + 1 */
                    if (0 == SearchFluidNodes(k, j, i, type, space)) {
                        space->nodeFlag[idx] = -space->nodeFlag[idx] + (type - 1) * geometry->totalN; /* geometry conserved */
                        break; /* exit search loop once identified successfully */
                    }
                }
            }
        }
    }
    return 0;
}
/*
 * Search around current node and check whether current node has at
 * least one fluid node regarding to the specified coordinate range.
 */
static int SearchFluidNodes(const int k, const int j, const int i, 
        const int h, const Space *space)
{
    const int idxW = IndexMath(k, j, i - h, space);
    const int idxE = IndexMath(k, j, i + h, space);
    const int idxS = IndexMath(k, j - h, i, space);
    const int idxN = IndexMath(k, j + h, i, space);
    const int idxF = IndexMath(k - h, j, i, space);
    const int idxB = IndexMath(k + h, j, i, space);
    return (space->nodeFlag[idxW] * space->nodeFlag[idxE] * 
            space->nodeFlag[idxS] * space->nodeFlag[idxN] * 
            space->nodeFlag[idxF] * space->nodeFlag[idxB]);
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
    Real nVec[DIMS] = {0.0}; /* normal vector */
    Real taVec[DIMS] = {0.0}; /* tangential vector */
    Real tbVec[DIMS] = {0.0}; /* tangential vector */
    Real rhs[DIMS] = {0.0}; /* right hand side vector */
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
                        InverseDistanceWeighting(UoImage, &weightSum, imageZ, imageY, imageX, 
                                ComputeK(imageZ, space), ComputeJ(imageY, space), ComputeI(imageX, space), 
                                1, FLUID, U, space, model, geometry);
                        /* enforce boundary conditions at boundary point */
                        if (0 < geo[GROUGH]) { /* noslip wall */
                            UoBC[1] = geo[GU];
                            UoBC[2] = geo[GV];
                            UoBC[3] = geo[GW];
                        } else { /* slip wall */
                            OrthogonalSpace(nVec, taVec, tbVec, info);
                            rhs[X] = geo[GU] * nVec[X] + geo[GV] * nVec[Y] + geo[GW] * nVec[Z];
                            rhs[Y] = (UoImage[1] * taVec[X] + UoImage[2] * taVec[Y] + UoImage[3] * taVec[Z]) / weightSum;
                            rhs[Z] = (UoImage[1] * tbVec[X] + UoImage[2] * tbVec[Y] + UoImage[3] * tbVec[Z]) / weightSum;
                            UoBC[1] = nVec[X] * rhs[X] + taVec[X] * rhs[Y] + tbVec[X] * rhs[Z];
                            UoBC[2] = nVec[Y] * rhs[X] + taVec[Y] * rhs[Y] + tbVec[Y] * rhs[Z];
                            UoBC[3] = nVec[Z] * rhs[X] + taVec[Z] * rhs[Y] + tbVec[Z] * rhs[Z];
                        }
                        UoBC[4] = UoImage[4] / weightSum; /* zero gradient pressure */
                        if (0 > geo[GT]) { /* adiabatic, dT/dn = 0 */
                            UoBC[5] = UoImage[5] / weightSum;
                        } else { /* otherwise, use specified constant wall temperature, T = Tw */
                            UoBC[5] = geo[GT];
                        }
                        /* add the boundary point as a stencil */
                        ApplyWeighting(UoImage, &weightSum, info[GSDS] * info[GSDS], UoBC, space->tinyL);
                        /* Normalize the weighted values of the image point */
                        NormalizeReconstructedValues(UoImage, weightSum);
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
                        /* Normalize the weighted values as reconstructed values. */
                        NormalizeReconstructedValues(UoGhost, weightSum);
                    }
                    UoGhost[0] = UoGhost[4] / (UoGhost[5] * model->gasR); /* compute density */
                    ConservativeByPrimitive(U, idx * DIMU, UoGhost, model);
                }
            }
        }
    }
    return 0;
}
int InverseDistanceWeighting(Real Uo[], Real *weightSum, const Real z, const Real y, const Real x,
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
                if (!((0 < Uoh[0]) && (FLT_MAX > Uoh[0]))) {
                    fprintf(stderr, "k=%d, j=%d, i=%d, rho=%.6g\n", k + kh, j + jh, i + ih, Uoh[0]);
                    FatalError("illegal density encountered, solution diverges...");
                }
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
int NormalizeReconstructedValues(Real Uo[], const Real weightSum)
{
    for (int n = 0; n < DIMUo; ++n) {
        Uo[n] = Uo[n] / weightSum;
    }
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
static int OrthogonalSpace(Real *nVec, Real *taVec, Real *tbVec, const Real *info)
{
    nVec[X] = info[GSNX]; 
    nVec[Y] = info[GSNY]; 
    nVec[Z] = info[GSNZ]; 
    int mark = Z; /* default mark for minimum component */
    if (fabs(nVec[mark]) > fabs(nVec[Y])) {
        mark = Y;
    }
    if (fabs(nVec[mark]) > fabs(nVec[X])) {
        mark = X;
    }
    if (X == mark) {
        taVec[X] = 0;
        taVec[Y] = -nVec[Z] / sqrt(nVec[Y] * nVec[Y] + nVec[Z] * nVec[Z]);
        taVec[Z] = nVec[Y] / sqrt(nVec[Y] * nVec[Y] + nVec[Z] * nVec[Z]);
    } else {
        if (Y == mark) {
            taVec[X] = -nVec[Z] / sqrt(nVec[X] * nVec[X] + nVec[Z] * nVec[Z]);
            taVec[Y] = 0;
            taVec[Z] = nVec[X] / sqrt(nVec[X] * nVec[X] + nVec[Z] * nVec[Z]);
        } else {
            taVec[X] = -nVec[Y] / sqrt(nVec[X] * nVec[X] + nVec[Y] * nVec[Y]);
            taVec[Y] = nVec[X] / sqrt(nVec[X] * nVec[X] + nVec[Y] * nVec[Y]);
            taVec[Z] = 0;
        }
    }
    tbVec[X] = taVec[Y] * nVec[Z] - nVec[Y] * taVec[Z];
    tbVec[Y] = -taVec[X] * nVec[Z] + nVec[X] * taVec[Z];
    tbVec[Z] = taVec[X] * nVec[Y] - nVec[X] * taVec[Y];
    return 0;
}
/* a good practice: end file with a newline */

