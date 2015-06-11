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
#include "solid_dynamics.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "gcibm.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int SurfaceForceIntegration(const Real *U, const Space *space,
        Particle *particle, const Partition *part, const Flow *flow);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ParticleSpatialEvolution(Real *U, const Real dt, Space *space, 
        Particle *particle, const Partition *part, const Flow *flow)
{
    /*
     * Compute the forces acting on particles.
     */
    SurfaceForceIntegration(U, space, particle, part, flow);
    /*
     * Update particle velocity and position
     */
    Real *ptk = NULL;
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = IndexGeometry(geoCount, particle);
        /* velocity: v(t) = v(t0) + f/m * dt */
        ptk[5] = ptk[5] + ptk[8] * ptk[13] * dt;
        ptk[6] = ptk[6] + ptk[9] * ptk[13] * dt;
        ptk[7] = ptk[7] + ptk[10] * ptk[13] * dt;
        /* spatial position: x(t) = x(t0) + v(t) * dt - 1/2 * f/m * dt^2 */
        ptk[0] = ptk[0] + ptk[5] * dt - 0.5 * ptk[8] * ptk[13] * dt * dt;
        ptk[1] = ptk[1] + ptk[6] * dt - 0.5 * ptk[9] * ptk[13] * dt * dt;
        ptk[2] = ptk[2] + ptk[7] * dt - 0.5 * ptk[10] * ptk[13] * dt * dt;
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
    int geoID = 0; /* geometry id */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (OFFSET > space->nodeFlag[idx]) { /* it's not a ghost */
                    continue;
                }
                geoID = space->nodeFlag[idx] - OFFSET; /* extract geometry information */
                ptk = IndexGeometry(geoID, particle);
                if (0 > InGeometry(k, j, i, ptk, space)) { /* still in the solid geometry */
                    continue;
                }
                /* reconstruction of flow values */
                InverseDistanceWeighting(Uo, ComputeZ(k, space), ComputeY(j, space), ComputeX(i, space), 
                        k, j, i, 2, U, space, flow);
                /* Normalize the weighted values as reconstructed values. */
                NormalizeReconstructedValues(Uo);
                ConservativeByPrimitive(U, idx * DIMU, Uo, flow);
            }
        }
    }
    /*
     * Recompute the changed geometry and apply boundary condition.
     */
    ComputeDomainGeometryGCIBM(space, particle, part);
    BoundaryConditionGCIBM(U, space, particle, part, flow);
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
        Particle *particle, const Partition *part, const Flow *flow)
{
    int idx = 0; /* linear array index math variable */
    int geoID = 0; /* geometry id */
    Real *ptk = NULL;
    Real info[INFOGEO] = {0.0}; /* store calculated geometry information */
    Real p = 0.0;
    /* reset some non accumulative information of particles to zero */
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = IndexGeometry(geoCount, particle);
        ptk[8] = 0; /* fx */
        ptk[9] = 0; /* fy */
        ptk[10] = 0; /* fz */
        ptk[11] = 0; /* tally */
    }
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (OFFSET > space->nodeFlag[idx]) { /* it's not a ghost */
                    continue;
                }
                geoID = space->nodeFlag[idx] - OFFSET; /* extract geometry information */
                ptk = IndexGeometry(geoID, particle);
                CalculateGeometryInformation(info, k, j, i, ptk, space);
                p = ComputePressure(idx * DIMU, U, flow);
                ptk[8] = ptk[8] - p * info[5]; /* integrate fx */
                ptk[9] = ptk[9] - p * info[6]; /* integrate fy */
                ptk[10] = ptk[10] - p * info[7]; /* integrate fz */
                ptk[11] = ptk[11] + 1; /* count number of ghosts */
            }
        }
    }
    /* calibrate the sum of discrete forces into integration */
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = IndexGeometry(geoCount, particle);
        ptk[8] = ptk[8] * ptk[12] / ptk[11];
        ptk[9] = ptk[9] * ptk[12] / ptk[11];
        ptk[10] = ptk[10] * ptk[12] / ptk[11];
    }
    return 0;
}
/* a good practice: end file with a newline */

