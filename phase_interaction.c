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
#include "phase_interaction.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "immersed_boundary.h"
#include "computational_geometry.h"
#include "linear_system.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void SurfaceForceIntegration(Space *, const Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int PhaseInteraction(const Real dt, Space *space, const Model *model)
{
    Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    /*
     * Compute the forces acting on particles.
     */
    SurfaceForceIntegration(space, model);
    /*
     * Update particle velocity and position
     */
    RealVec a = {0.0}; /* acceleration */
    RealVec alpha = {0.0}; /* angular acceleration */
    RealVec offset = {0.0}; /* translation */
    RealVec angle = {0.0}; /* rotation */
    const RealVec scale = {1.0, 1.0, 1.0}; /* scale */
    for (int n = 0; n < geo->totalN; ++n) {
        poly = geo->poly + n;
        MatrixLinearSystemSolver(DIMS, poly->I, 1, (Real (*)[1])alpha, (Real (*)[1])poly->Tau);
        for (int s = 0; s < DIMS; ++s) {
            /* acceleration from surface force and body force */
            a[s] = poly->F[s] / (poly->rho * poly->volume) + poly->gState * model->g[s];
            alpha[s] = alpha[s] / poly->rho;
        }
        for (int s = 0; s < DIMS; ++s) {
            /* velocity integration v(t[n+1]) = v(t[n]) + a * dt */
            poly->V[s] = poly->V[s] + a[s] * dt;
            poly->W[s] = poly->W[s] + alpha[s] * dt;
            /* position integration x(t[n+1]) = x(t[n]) + v(t[n+1]) * dt - 1/2 * a * dt^2 */
            offset[s] = poly->V[s] * dt - 0.5 * a[s] * dt * dt;
            angle[s] = poly->W[s] * dt - 0.5 * alpha[s] * dt * dt;
        }
        if (0 == poly->faceN) { /* analytical sphere */
            poly->O[X] = poly->O[X] + offset[X];
            poly->O[Y] = poly->O[Y] + offset[Y];
            poly->O[Z] = poly->O[Z] + offset[Z];
        } else { /* triangulated polyhedron */
            Transformation(scale, angle, offset, poly);
        }
    }
    /*
     * Recompute the changed geometry and apply boundary condition.
     */
    ComputeGeometryDomain(space, model);
    ImmersedBoundaryTreatment(TO, space, model);
    return 0;
}
static void SurfaceForceIntegration(Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    const Node *node = space->node;
    const Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    Real p = 0.0;
    /* reset some non accumulative information of particles to zero */
    for (int n = 0; n < geo->totalN; ++n) {
        poly = geo->poly + n;
        poly->F[X] = 0; /* fx */
        poly->F[Y] = 0; /* fy */
        poly->F[Z] = 0; /* fz */
        poly->Tau[X] = 0; /* torque x */
        poly->Tau[Y] = 0; /* torque y */
        poly->Tau[Z] = 0; /* torque z */
        poly->nodeN = 0; /* tally */
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

