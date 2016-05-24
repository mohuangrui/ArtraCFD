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
#include "solve.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <float.h> /* size of floating point values */
#include "initialization.h"
#include "fluid_dynamics.h"
#include "phase_interaction.h"
#include "data_stream.h"
#include "timer.h"
#include "data_probe.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int SolutionEvolution(Time *, Space *, const Model *);
static Real ComputeTimeStep(const Time *, const Space *, const Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int Solve(Time *time, Space *space, const Model *model)
{
    ShowInformation("Solving...");
    fprintf(stdout, "  initializing...\n");
    InitializeComputationalDomain(time, space, model);
    fprintf(stdout, "  time marching...\n");
    SolutionEvolution(time, space, model);
    ShowInformation("Session End");
    return 0;
}
static int SolutionEvolution(Time *time, Space *space, const Model *model)
{
    Real dt = time->end - time->now; /* time step size */
    /* check whether current time is equal to or larger than the end time */
    if (0.0 >= dt) {
        fprintf(stdout, "  current time is equal to or larger than the end time...\n");
        return 1;
    }
    /* obtain the desired export time interval */
    const Real interval = dt / (Real)(MaxInt(time->outputN - time->countOutput, 1));
    const Real probeInterval = dt / (Real)(time->outputNProbe);
    Real record = 0.0; /* used for control when to export data */
    Real probeRecord = 0.0; /* used for control when to export probe data */
    /* set some timers for monitoring time consuming of process */
    Timer operationTimer; /* timer for computing operations */
    Real operationTime = 0.0; /* record consuming time of operation */
    while ((time->now < time->end) && (time->countStep <= time->stepN)) {
        /*
         * Step count
         */
        ++(time->countStep);
        /*
         * Calculate dt for current time step
         */
        dt = ComputeTimeStep(time, space, model);
        /*
         * Update current time stamp, if current time exceeds the end time, 
         * recompute the value of dt to make current time equal to the end time.
         */
        time->now = time->now + dt;
        if (time->now > time->end) {
            dt = time->end - (time->now - dt);
            time->now = time->end;
        }
        fprintf(stdout, "\nstep=%d; time=%.6g; remain=%.6g; dt=%.6g;\n", 
                time->countStep, time->now, time->end - time->now, dt);
        /*
         * Compute field data in current time step, treat phase interaction
         * as two physical processes, and split these two processes in time
         * space use a technique similar to Strang's splitting method.
         */
        TickTime(&operationTimer);
        if (1 == model->fsi) {
            PhaseInteraction(0.5 * dt, space, model);
        }
        FluidDynamics(dt, space, model);
        if (1 == model->fsi) {
            PhaseInteraction(0.5 * dt, space, model);
        }
        operationTime = TockTime(&operationTimer);
        fprintf(stdout, "  elapsed: %.6gs\n", operationTime);
        /*
         * Export data if accumulated time increases to anticipated interval.
         */
        record = record + dt;
        probeRecord = probeRecord + dt;
        if ((record >= interval) || (time->now >= time->end) || (time->countStep == time->stepN)) {
            ++(time->countOutput); /* export count increase */
            fprintf(stdout, "  exporting data...\n");
            WriteFieldData(time, space, model);
            WriteGeometryData(time, &(space->geo));
            record = 0.0; /* reset accumulated time */
        }
        if ((probeRecord > probeInterval) || (time->now >= time->end) || (time->countStep == time->stepN)) {
            WriteFieldDataAtProbes(time, space, model);
            probeRecord = 0.0; /* reset probe accumulated time */
        }
    }
    return 0;
}
static Real ComputeTimeStep(const Time *time, const Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    const Geometry *geo = &(space->geo);
    const Polyhedron *poly = NULL;
    const Node *const node = space->node;
    const Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    Real speed = 0.0;
    Real speedMax = FLT_MIN;
    /*
     * Incorporate solid dynamics into CFL condition.
     */
    for (int n = 0; n < geo->totalN; ++n) {
        poly = geo->poly + n;
        speed = MaxReal(fabs(poly->V[X]) + fabs(poly->W[X]) * poly->r, 
                MaxReal(fabs(poly->V[Y]) + fabs(poly->W[Y]) * poly->r, 
                    fabs(poly->V[Z]) + fabs(poly->W[Z]) * poly->r));
        if (speedMax < speed) {
            speedMax = speed;
        }
    }
    /*
     * Incorporate fluid dynamics into CFL condition.
     */
    Real Uo[DIMUo] = {0.0};
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                U = node[idx].U[TO];
                if (0 != node[idx].geoID) {
                    continue;
                }
                PrimitiveByConservative(model->gamma, model->gasR, U, Uo);
                speed = MaxReal(fabs(Uo[1]), MaxReal(fabs(Uo[2]), fabs(Uo[3]))) + 
                    sqrt(model->gamma * model->gasR * Uo[5]);
                if (speedMax < speed) {
                    speedMax = speed;
                }
            }
        }
    }
    return time->numCFL * MinReal(part->d[X], MinReal(part->d[Y], part->d[Z])) / speedMax;
}
/* a good practice: end file with a newline */

