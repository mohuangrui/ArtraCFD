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
#include "fluid_solid_interaction.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "gcibm.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int SurfaceForceIntegration(const Real *, const Space *,
        const Model *, const Partition *, Geometry *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int FluidSolidInteraction(Real *U, Space *space, const Model *model, const Partition *part,
        Geometry *geometry, const Real dt)
{
    /*
     * Compute the forces acting on particles.
     */
    SurfaceForceIntegration(U, space, model, part, geometry);
    /*
     * Update particle velocity and position
     */
    Real *geo = NULL;
    for (int geoCount = 0; geoCount < geometry->totalN; ++geoCount) {
        geo = IndexGeometry(geoCount, geometry);
        /* velocity: v(t[n+1]) = v(t[n]) + f * (1/m) * dt */
        geo[GU] = geo[GU] + geo[GFX] * (1/geo[GMASS]) * dt;
        geo[GV] = geo[GV] + geo[GFY] * (1/geo[GMASS]) * dt;
        geo[GW] = geo[GW] + geo[GFZ] * (1/geo[GMASS]) * dt;
        /* spatial position: x(t[n+1]) = x(t[n]) + v(t[n+1]) * dt - 1/2 * f * (1/m) * dt^2 */
        geo[GX] = geo[GX] + geo[GU] * dt - 0.5 * geo[GFX] * (1/geo[GMASS]) * dt * dt;
        geo[GY] = geo[GY] + geo[GV] * dt - 0.5 * geo[GFY] * (1/geo[GMASS]) * dt * dt;
        geo[GZ] = geo[GZ] + geo[GW] * dt - 0.5 * geo[GFZ] * (1/geo[GMASS]) * dt * dt;
    }
    /*
     * After the spatial positions of particles updated, some inner nodes fall
     * out of the solid region, these nodes need to be identified and the flow
     * values of these nodes need to be interpolated. Since the speed of
     * particles will not greater than the speed of fluid, the spatial motions
     * of particles will not exceed one grid per time step based on the CFL
     * condition of flow. Therefore, only ghost nodes need to be verified.
     */
    int idx = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0}; /* store weighted primitives */
    Real UoBC[DIMUo] = {0.0}; /* physical primitives at boundary point */
    Real info[INFOGHOST] = {0.0}; /* store calculated geometry information */
    int geoID = 0; /* geometry id */
    int type = 0;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (OFFSET > space->nodeFlag[idx]) { /* it's not a ghost */
                    continue;
                }
                for (type = 1; type < space->ng + 2; ++type) { /* extract geometry identifier */
                    geoID = space->nodeFlag[idx] - OFFSET - (type - 1) * geometry->totalN;
                    if ((0 <= geoID) && (geometry->totalN > geoID)) { /* a ghost node with current type */
                        break;
                    }
                }
                geo = IndexGeometry(geoID, geometry);
                if (0 > InGeometry(k, j, i, geo, space)) { /* still in the solid geometry */
                    continue;
                }
                /* reconstruction of flow values */
                CalculateGeometryInformation(info, k, j, i, geo, space);
                FlowReconstruction(Uo, info[GSZ], info[GSY], info[GSX], k, j, i, type,
                        UoBC, info, geo, U, space, model, geometry);
                Uo[0] = Uo[4] / (Uo[5] * model->gasR); /* compute density */
                ConservativeByPrimitive(U, idx * DIMU, Uo, model);
            }
        }
    }
    /*
     * Recompute the changed geometry and apply boundary condition.
     */
    ComputeGeometryDomain(space, part, geometry);
    BoundaryTreatmentsGCIBM(U, space, model, part, geometry);
    return 0;
}
/*
 * Add the pressure force at boundary to the pressure force of corresponding
 * particle. The pressure at boundary will equal to the pressure of the ghost
 * node since zero pressure gradient at wall normal direction is enforced here.
 * A even spaced pressure distribution over the particle surface is assumed 
 * since we only compute the pressure at the boundary point that has a ghost 
 * neighbour. By this approach, the accuracy of pressure integration along 
 * particle surface will increase correspondingly with the increase of mesh 
 * resolution while saving remarkable computation effort.
 */
static int SurfaceForceIntegration(const Real *U, const Space *space,
        const Model *model, const Partition *part, Geometry *geometry)
{
    int idx = 0; /* linear array index math variable */
    int geoID = 0; /* geometry id */
    Real *geo = NULL;
    Real info[INFOGHOST] = {0.0}; /* store calculated geometry information */
    Real p = 0.0;
    /* reset some non accumulative information of particles to zero */
    for (int geoCount = 0; geoCount < geometry->totalN; ++geoCount) {
        geo = IndexGeometry(geoCount, geometry);
        geo[GFX] = 0; /* fx */
        geo[GFY] = 0; /* fy */
        geo[GFZ] = 0; /* fz */
        geo[GTALLY] = 0; /* tally */
    }
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (OFFSET > space->nodeFlag[idx]) { /* it's not a ghost */
                    continue;
                }
                geoID = space->nodeFlag[idx] - OFFSET;
                if ((0 > geoID) || (geometry->totalN <= geoID)) { /* not a first type ghost node */
                    continue;
                }
                geo = IndexGeometry(geoID, geometry);
                CalculateGeometryInformation(info, k, j, i, geo, space);
                p = ComputePressure(idx * DIMU, U, model);
                geo[GFX] = geo[GFX] - p * info[GSNX]; /* integrate fx */
                geo[GFY] = geo[GFY] - p * info[GSNY]; /* integrate fy */
                geo[GFZ] = geo[GFZ] - p * info[GSNZ]; /* integrate fz */
                geo[GTALLY] = geo[GTALLY] + 1; /* count number of ghosts */
            }
        }
    }
    /* calibrate the sum of discrete forces into integration */
    for (int geoCount = 0; geoCount < geometry->totalN; ++geoCount) {
        geo = IndexGeometry(geoCount, geometry);
        geo[GFX] = geo[GFX] * geo[GAREA] / geo[GTALLY];
        geo[GFY] = geo[GFY] * geo[GAREA] / geo[GTALLY];
        geo[GFZ] = geo[GFZ] * geo[GAREA] / geo[GTALLY];
    }
    return 0;
}
/* a good practice: end file with a newline */

