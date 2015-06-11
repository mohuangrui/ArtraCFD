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
#include "solve.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <float.h> /* size of floating point values */
#include "initialization.h"
#include "temporal_discretization.h"
#include "particle_dynamics.h"
#include "data_stream.h"
#include "geometry_stream.h"
#include "timer.h"
#include "probe.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int SolutionEvolution(Field *, Space *, Particle *, Time *, 
        const Partition *, const Flow *);
static Real ComputeTimeStep(const Real *U, const Space *, 
        const Particle *, const Time *, const Partition *, const Flow *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int Solve(Field *field, Space *space, Particle *particle, Time *time, 
        const Partition *part, const Flow *flow)
{
    InitializeFlowField(field->U, space, particle, time, part, flow);
    SolutionEvolution(field, space, particle, time, part, flow);
    return 0;
}
static int SolutionEvolution(Field *field, Space *space, Particle *particle,
        Time *time, const Partition *part, const Flow *flow)
{
    ShowInformation("Time marching...");
    /* check whether current time is equal to or larger than the total time */
    if (0 >= (time->totalTime - time->currentTime)) {
        ShowInformation("  current time is equal to or larger than total time...");
        ShowInformation("Session End");
        return 1;
    }
    /* obtain the desired export time interval */
    const Real exportTimeInterval = (time->totalTime - time->currentTime) / 
        (Real)(MaxInt(time->totalOutputTimes - time->outputCount, 1));
    const Real probeExportInterval = (time->totalTime - time->currentTime) / 
        (Real)(flow->outputProbe);
    Real accumulatedTime = 0.0; /* used for control when to export data */
    Real probeAccumulatedTime = 0.0; /* used for control when to export probe data */
    /* set some timers for monitoring time consuming of process */
    Timer operationTimer; /* timer for computing operations */
    Real operationTime = 0.0; /* record consuming time of operation */
    for (time->stepCount += 1; (time->currentTime < time->totalTime) && 
            (time->stepCount <= time->totalStep); ++(time->stepCount)) {
        /*
         * Calculate dt for current time step
         */
        time->dt = ComputeTimeStep(field->U, space, particle, time, part, flow);
        /*
         * Update current time stamp, if current time exceeds the total time, 
         * recompute the value of dt to make current time equal total time.
         */
        time->currentTime = time->currentTime + time->dt;
        if (time->currentTime > time->totalTime) { /* need to refine "dt" to reach totTime  */
            time->dt = time->totalTime - (time->currentTime - time->dt);
            time->currentTime = time->totalTime;
        }
        fprintf(stdout, "\nstep=%d; time=%.6g; remain=%.6g; dt=%.6g;\n", time->stepCount, 
                time->currentTime, time->totalTime - time->currentTime, time->dt);
        /*
         * Compute field data in current time step, treat the interaction of
         * fluid and particles as two physical processes, and split these two
         * processes in time space use a technique similar to Strang's
         * splitting method.
         */
        TickTime(&operationTimer);
        ParticleDynamics(field->U, 0.5 * time->dt, space, particle, part, flow);
        FluidDynamics(field, time->dt, space, particle, part, flow);
        ParticleDynamics(field->U, 0.5 * time->dt, space, particle, part, flow);
        operationTime = TockTime(&operationTimer);
        fprintf(stdout, "  elapsed: %.6gs\n", operationTime);
        /*
         * Export computed data. Use accumulatedTime as a flag, if
         * accumulatedTime increases to anticipated export interval,
         * write data out. Because the accumulatedTime very likely
         * can not increase to the exporting interval at the last
         * phase, then add a extra condition that if current time
         * is the total time, then also write data out.
         */
        accumulatedTime = accumulatedTime + time->dt;
        probeAccumulatedTime = probeAccumulatedTime + time->dt;
        if ((accumulatedTime >=  exportTimeInterval) || 
                (0 == time->currentTime - time->totalTime) ||
                (time->stepCount == time->totalStep)) {
            ++(time->outputCount); /* export count increase */
            TickTime(&operationTimer);
            WriteComputedData(field->U, space, time, part, flow);
            WriteGeometryData(particle, time);
            operationTime = TockTime(&operationTimer);
            accumulatedTime = 0; /* reset accumulated time */
            fprintf(stdout, "  data export time consuming: %.6gs\n", operationTime);
        }
        if ((probeAccumulatedTime >=  probeExportInterval) || 
                (0 == time->currentTime - time->totalTime) ||
                (time->stepCount == time->totalStep)) {
            WriteComputedDataAtProbes(time->stepCount, field->U, space, part, flow);
            probeAccumulatedTime = 0; /* reset probe accumulated time */
        }
    }
    ShowInformation("Session End");
    return 0;
}
static Real ComputeTimeStepByCFL(const Real *U, const Space *space, const Particle *particle, 
        const Time *time, const Partition *part, const Flow *flow)
{
    Real velocity = 0.0;
    Real velocityMax = FLT_MIN;
    /*
     * Incorporate particle dynamics into CFL condition.
     */
    Real *ptk = NULL;
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = IndexGeometry(geoCount, particle);
        velocity = MaxReal(fabs(ptk[7]), MaxReal(fabs(ptk[5]), fabs(ptk[6])));
        if (velocityMax < velocity) {
            velocityMax = velocity;
        }
    }
    /*
     * Incorporate fluid dynamics into CFL condition.
     */
    int idx = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0};
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (FLUID != space->nodeFlag[idx]) {
                    continue;
                }
                PrimitiveByConservative(Uo, idx * DIMU, U, flow);
                velocity = MaxReal(fabs(Uo[3]), MaxReal(fabs(Uo[1]), fabs(Uo[2]))) + sqrt((Uo[4] / Uo[0]) * flow->gamma);
                if (velocityMax < velocity) {
                    velocityMax = velocity;
                }
            }
        }
    }
    return time->numCFL * MinReal(space->dx, MinReal(space->dz, space->dy)) / velocityMax;
}
/* a good practice: end file with a newline */

