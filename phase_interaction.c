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
    const Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    RealVec pG = {0.0}; /* ghost point */
    RealVec pO = {0.0}; /* boundary point */
    RealVec pI = {0.0}; /* image point */
    RealVec N = {0.0}; /* normal */
    Real p = {0.0};
    RealVec r = {0.0}; /* position vector */
    RealVec F = {0.0}; /* force */
    RealVec Tau = {0.0}; /* touque */
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
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                if (0 == node[idx].geoID) {
                    continue;
                }
                if (1 != node[idx].ghostID) {
                    continue;
                }
                poly = geo->poly + node[idx].geoID - 1;
                ComputeGeometricData(node[idx].faceID, poly, pG, pO, pI, N);
                p = ComputePressure(model->gamma, node[idx].U[TO]);
                r[X] = pO[X] - poly->O[X];
                r[Y] = pO[Y] - poly->O[Y];
                r[Z] = pO[Z] - poly->O[Z];
                F[X] = -p * N[X];
                F[Y] = -p * N[Y];
                F[Z] = -p * N[Z];
                Cross(r, F, Tau);
                poly->F[X] = poly->F[X] + F[X]; /* integrate fx */
                poly->F[Y] = poly->F[Y] + F[Y]; /* integrate fy */
                poly->F[Z] = poly->F[Z] + F[Z]; /* integrate fz */
                poly->Tau[X] = poly->Tau[X] + Tau[X]; /* integrate torque x */
                poly->Tau[Y] = poly->Tau[Y] + Tau[Y]; /* integrate torque y */
                poly->Tau[Z] = poly->Tau[Z] + Tau[Z]; /* integrate torque z */
                ++(poly->nodeN); /* count number of ghosts */
            }
        }
    }
    /* calibrate the sum of discrete forces into integration */
    for (int n = 0; n < geo->totalN; ++n) {
        poly = geo->poly + n;
        poly->F[X] = poly->F[X] * poly->area / poly->nodeN;
        poly->F[Y] = poly->F[Y] * poly->area / poly->nodeN;
        poly->F[Z] = poly->F[Z] * poly->area / poly->nodeN;
        poly->Tau[X] = poly->Tau[X] * poly->area / poly->nodeN;
        poly->Tau[Y] = poly->Tau[Y] * poly->area / poly->nodeN;
        poly->Tau[Z] = poly->Tau[Z] * poly->area / poly->nodeN;
    }
    return;
}
/* a good practice: end file with a newline */

