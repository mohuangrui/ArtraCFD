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
#include "solid_dynamics.h"
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
/*
 * Mo, H., Lien, F. S., Zhang, F., & Cronin, D. S. (2017). A numerical
 * framework for the direct simulation of dense particulate flow under
 * explosive dispersal. Shock Waves, 1-19.
 */
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
    Real dt = time->end - time->now;
    const Real zero = 0.0;
    /* check whether current time is equal to or larger than the end time */
    if (zero >= dt) {
        fprintf(stdout, "  current time is equal to or larger than the end time...\n");
        return 1;
    }
    /* obtain the desired data writing interval */
    const Real dtField = time->end / (Real)(time->writeN);
    const Real dtPoint = time->end / (Real)(time->pointWriteN);
    const Real dtLine = time->end / (Real)(time->lineWriteN);
    const Real dtCurve = time->end / (Real)(time->curveWriteN);
    const Real dtForce = time->end / (Real)(time->forceWriteN);
    Real recField = zero; /* field data writing recorder */
    Real recPoint = zero; /* point probe data writing recorder */
    Real recLine = zero; /* line probe data writing recorder */
    Real recCurve = zero; /* curve probe data writing recorder */
    Real recForce = zero; /* force probe data writing recorder */
    /* set some timers for monitoring time consuming of process */
    Timer timer; /* timer for computing operations */
    while ((time->now < time->end) && (time->stepC < time->stepN)) {
        /*
         * Step count
         */
        ++(time->stepC);
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
                time->stepC, time->now, time->end - time->now, dt);
        /*
         * Compute field data in current time step, treat phase interaction
         * as two physical processes, and split these two processes in time
         * space use Strang's splitting method.
         */
        TickTime(&timer);
        if (0 != model->fsi) {
            SolidDynamics(time->now, 0.5 * dt, space, model);
        }
        FluidDynamics(dt, space, model);
        if (0 != model->fsi) {
            SolidDynamics(time->now, 0.5 * dt, space, model);
        }
        fprintf(stdout, "  elapsed: %.6gs\n", TockTime(&timer));
        /*
         * Export data if accumulated time increases to anticipated interval.
         */
        recField = recField + dt;
        recPoint = recPoint + dt;
        recLine = recLine + dt;
        recCurve = recCurve + dt;
        recForce = recForce + dt;
        if ((recForce > dtForce) || (time->now == time->end) || (time->stepC == time->stepN)) {
            SurfaceForceIntegration(space, model);
            WriteSurfaceForceData(time, space);
            recForce = zero; /* reset probe accumulated time */
        }
        if ((recField > dtField) || (time->now == time->end) || (time->stepC == time->stepN)) {
            ++(time->writeC); /* export count increase */
            fprintf(stdout, "  writing field data...\n");
            WriteFieldData(time, space, model);
            WriteGeometryData(time, &(space->geo));
            recField = zero; /* reset accumulated time */
        }
        if ((recPoint > dtPoint) || (time->now == time->end) || (time->stepC == time->stepN)) {
            WriteFieldDataAtPointProbes(time, space, model);
            recPoint = zero; /* reset probe accumulated time */
        }
        if ((recLine > dtLine) || (time->now == time->end) || (time->stepC == time->stepN)) {
            WriteFieldDataAtLineProbes(time, space, model);
            recLine = zero; /* reset probe accumulated time */
        }
        if ((recCurve > dtCurve) || (time->now == time->end) || (time->stepC == time->stepN)) {
            WriteFieldDataAtCurveProbes(time, space, model);
            recCurve = zero; /* reset probe accumulated time */
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
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        speed = MaxReal(fabs(poly->V[TO][X]), MaxReal(fabs(poly->V[TO][Y]), fabs(poly->V[TO][Z]))) +
            MaxReal(fabs(poly->W[TO][X]), MaxReal(fabs(poly->W[TO][Y]), fabs(poly->W[TO][Z]))) * poly->r;
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
                if (0 != node[idx].gid) {
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

