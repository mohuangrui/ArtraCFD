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
#include <limits.h> /* sizes of integral types */
#include "initialization.h"
#include "fluid_dynamics.h"
#include "solid_dynamics.h"
#include "data_stream.h"
#include "timer.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void EvolveSolution(Time *, Space *, const Model *);
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
    ShowInfo("Solving...\n");
    ShowInfo("  initializing...\n");
    InitializeComputeDomain(time, space, model);
    ShowInfo("  time marching...\n");
    EvolveSolution(time, space, model);
    ShowInfo("Session");
    return 0;
}
static void EvolveSolution(Time *time, Space *space, const Model *model)
{
    Real dt = time->end - time->now;
    const Real zero = 0.0;
    if (zero >= dt) {
        ShowWarning("  time.now >= time.end");
        return;
    }
    Timer tm; /* timer for computing operations */
    /* data writing interval and recorder */
    const Real dtData[NPROBE] = {time->end / (Real)(time->dataW[PROPT]),
        time->end / (Real)(time->dataW[PROLN]), time->end / (Real)(time->dataW[PROCV]),
        time->end / (Real)(time->dataW[PROFC]), time->end / (Real)(time->dataW[PROSD])};
    Real rcData[NPROBE] = {zero};
    /* time instants interval and recorder */
    const Real tmInt = (INT_MAX == time->dataW[PROSD]) ? time->end : dtData[PROSD]; /* a specific instant */
    Real rcInt = zero; /* time instant recorder */
    while ((time->now < time->end) && (time->stepC < time->stepN)) {
        ++(time->stepC);
        dt = ComputeTimeStep(time, space, model);
        if (rcInt + dt > tmInt) { /* rectify dt */
            dt = tmInt - rcInt;
            rcInt = zero;
        } else {
            rcInt = rcInt + dt;
        }
        time->now = time->now + dt;
        if (time->now > time->end) { /* rectify dt */
            dt = time->end - (time->now - dt);
            time->now = time->end;
        }
        ShowInfo("\nstep=%d; time=%.6g; remain=%.6g; dt=%.6g;\n",
                time->stepC, time->now, time->end - time->now, dt);
        TickTime(&tm);
        if (0 != model->psi) {
            EvolveSolidDynamics(time->now, 0.5 * dt, space, model);
        }
        EvolveFluidDynamics(dt, space, model);
        if (0 != model->psi) {
            EvolveSolidDynamics(time->now, 0.5 * dt, space, model);
        }
        ShowInfo("  elapsed: %.6gs\n", TockTime(&tm));
        /* export data if accumulated time increases to anticipated interval */
        for (int n = 0; n < NPROBE; ++n) {
            rcData[n] = rcData[n] + dt;
            if ((rcData[n] >= dtData[n]) || (time->now == time->end) || (time->stepC == time->stepN)) {
                if (PROFC == n) {
                    IntegrateSurfaceForce(space, model);
                }
                if (PROSD == n) {
                    ShowInfo("  writing data...\n");
                    ++(time->dataC); /* export count increase */
                }
                WriteData(n, time, space, model);
                rcData[n] = zero; /* reset probe accumulated time */
            }
        }
    }
    return;
}
static Real ComputeTimeStep(const Time *time, const Space *space, const Model *model)
{
    const Partition *const part = &(space->part);
    const Node *const node = space->node;
    const Geometry *const geo = &(space->geo);
    const Polyhedron *poly = NULL;
    const Real *restrict U = NULL;
    Real Uo[DIMUo] = {0.0};
    int idx = 0; /* linear array index math variable */
    Real c = 0.0; /* speed of sound */
    RealVec V = {0.0}; /* characteristic speeds in each direction */
    RealVec Vmax = {0.0}; /* maximum characteristic speeds in each direction */
    /* incorporate solid dynamics into CFL condition */
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        V[X] = fabs(poly->V[TO][X]) + MaxReal(fabs(poly->W[TO][Y]), fabs(poly->W[TO][Z])) * poly->r;
        V[Y] = fabs(poly->V[TO][Y]) + MaxReal(fabs(poly->W[TO][Z]), fabs(poly->W[TO][X])) * poly->r;
        V[Z] = fabs(poly->V[TO][Z]) + MaxReal(fabs(poly->W[TO][X]), fabs(poly->W[TO][Y])) * poly->r;
        for (int s = 0; s < DIMS; ++s) {
            if (Vmax[s] < V[s]) {
                Vmax[s] = V[s];
            }
        }
    }
    /* incorporate fluid dynamics into CFL condition */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                U = node[idx].U[TO];
                if (0 != node[idx].did) {
                    continue;
                }
                MapPrimitive(model->gamma, model->gasR, U, Uo);
                c = sqrt(model->gamma * model->gasR * Uo[5]);
                for (int s = 0; s < DIMS; ++s) {
                    V[s] = fabs(Uo[s+1]) + c;
                    if (Vmax[s] < V[s]) {
                        Vmax[s] = V[s];
                    }
                }
            }
        }
    }
    return time->numCFL * MinReal(part->d[X] / Vmax[X], MinReal(part->d[Y] / Vmax[Y], part->d[Z] / Vmax[Z]));
}
/* a good practice: end file with a newline */

