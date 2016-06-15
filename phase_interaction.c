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
#include <float.h> /* size of floating point values */
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
int PhaseInteraction(const Real now, const Real dt, Space *space, const Model *model)
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
        if (1 == poly->state) { /* stationary object */
            continue;
        }
        if (now > poly->to) { /* end power supply */
            poly->a[X] = 0.0;
            poly->a[Y] = 0.0;
            poly->a[Z] = 0.0;
            poly->alpha[X] = 0.0;
            poly->alpha[Y] = 0.0;
            poly->alpha[Z] = 0.0;
            poly->to = FLT_MAX; /* avoid repeating */
        }
        /* tranlation and angular acceleration */
        MatrixLinearSystemSolver(DIMS, poly->I, 1, (Real (*)[1])alpha, (Real (*)[1])poly->Tau);
        for (int s = 0; s < DIMS; ++s) {
            /* acceleration from surface force and body force */
            a[s] = (poly->Fp[s] + poly->Fv[s]) / (poly->rho * poly->volume) + poly->a[s] + poly->g[s] + poly->ac[s];
            alpha[s] = alpha[s] / poly->rho + poly->alpha[s];
        }
        /* velocity and position integration */
        for (int s = 0; s < DIMS; ++s) {
            /* velocity integration v(t[n+1]) = v(t[n]) + a * dt */
            poly->V[s] = poly->V[s] + a[s] * dt;
            poly->W[s] = poly->W[s] + alpha[s] * dt;
            /* position integration x(t[n+1]) = x(t[n]) + v(t[n+1]) * dt - 1/2 * a * dt^2 */
            offset[s] = poly->V[s] * dt - 0.5 * a[s] * dt * dt;
            angle[s] = poly->W[s] * dt - 0.5 * alpha[s] * dt * dt;
        }
        /* transform geometry */
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
    const IntVec nMin = {part->ns[PIN][X][MIN], part->ns[PIN][Y][MIN], part->ns[PIN][Z][MIN]};
    const IntVec nMax = {part->ns[PIN][X][MAX], part->ns[PIN][Y][MAX], part->ns[PIN][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const int ng = part->ng;
    int idx = 0; /* linear array index math variable */
    int box[DIMS][LIMIT] = {{0}}; /* bounding box in node space */
    int nodeN = 0; /* count total number of interfacial nodes */
    const Real zero = 0.0;
    RealVec pG = {zero}; /* ghost point */
    RealVec pO = {zero}; /* boundary point */
    RealVec pI = {zero}; /* image point */
    RealVec N = {zero}; /* normal */
    Real Uo[DIMUo] = {zero};
    RealVec V = {zero};
    RealVec r = {zero}; /* position vector */
    RealVec Fp = {zero}; /* pressure force */
    RealVec Fv = {zero}; /* viscous force */
    RealVec Fs = {zero}; /* surface force */
    RealVec Tau = {zero}; /* torque */
    Real mu = zero; /* viscosity */
    for (int n = 0; n < geo->totalN; ++n) {
        poly = geo->poly + n;
        if (0 < poly->state) { /* surface force negligible */
            continue;
        }
        /* reset some non accumulative information to zero */
        poly->Fp[X] = zero;
        poly->Fp[Y] = zero;
        poly->Fp[Z] = zero;
        poly->Fv[X] = zero;
        poly->Fv[Y] = zero;
        poly->Fv[Z] = zero;
        poly->Tau[X] = zero;
        poly->Tau[Y] = zero;
        poly->Tau[Z] = zero;
        nodeN = 0; /* tally */
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ValidNodeSpace(NodeSpace(poly->box[s][MIN], sMin[s], dd[s], ng), nMin[s], nMax[s]);
            box[s][MAX] = ValidNodeSpace(NodeSpace(poly->box[s][MAX], sMin[s], dd[s], ng), nMin[s], nMax[s]) + 1;
        }
        for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
            for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    if ((2 == node[idx].layerID) && (n + 1 == node[idx].geoID)) {
                        ++nodeN; /* an interfacial node of current geometry */
                    }
                    if ((2 != node[idx].ghostID) || (n + 1 != node[idx].geoID)) {
                        continue;
                    }
                    /* surface force exerted by fluid (pressure + shear force) */
                    pG[X] = PointSpace(i, sMin[X], d[X], ng);
                    pG[Y] = PointSpace(j, sMin[Y], d[Y], ng);
                    pG[Z] = PointSpace(k, sMin[Z], d[Z], ng);
                    ComputeGeometricData(node[idx].faceID, poly, pG, pO, pI, N);
                    r[X] = pO[X] - poly->O[X];
                    r[Y] = pO[Y] - poly->O[Y];
                    r[Z] = pO[Z] - poly->O[Z];
                    PrimitiveByConservative(model->gamma, model->gasR, node[idx].U[TO], Uo);
                    Fp[X] = Uo[4] * N[X];
                    Fp[Y] = Uo[4] * N[Y];
                    Fp[Z] = Uo[4] * N[Z];
                    if ((zero < model->refMu) && (zero < poly->cf)) {
                        mu = model->refMu * Viscosity(Uo[5] * model->refT);
                        Cross(poly->W, r, V);
                        V[X] = poly->V[X] + V[X];
                        V[Y] = poly->V[Y] + V[Y];
                        V[Z] = poly->V[Z] + V[Z];
                        V[X] = Uo[1] - V[X];
                        V[Y] = Uo[2] - V[Y];
                        V[Z] = Uo[3] - V[Z];
                        Fv[X] = mu * (V[X] - Dot(V,N) * N[X]) / Dist(pG, pO);
                        Fv[Y] = mu * (V[Y] - Dot(V,N) * N[Y]) / Dist(pG, pO);
                        Fv[Z] = mu * (V[Z] - Dot(V,N) * N[Z]) / Dist(pG, pO);
                    } else {
                        Fv[X] = zero;
                        Fv[Y] = zero;
                        Fv[Z] = zero;
                    }
                    Fs[X] = Fp[X] + Fv[X];
                    Fs[Y] = Fp[Y] + Fv[Y];
                    Fs[Z] = Fp[Z] + Fv[Z];
                    Cross(r, Fs, Tau);
                    /* integration sum */
                    poly->Fp[X] = poly->Fp[X] + Fp[X];
                    poly->Fp[Y] = poly->Fp[Y] + Fp[Y];
                    poly->Fp[Z] = poly->Fp[Z] + Fp[Z];
                    poly->Fv[X] = poly->Fv[X] + Fv[X];
                    poly->Fv[Y] = poly->Fv[Y] + Fv[Y];
                    poly->Fv[Z] = poly->Fv[Z] + Fv[Z];
                    poly->Tau[X] = poly->Tau[X] + Tau[X];
                    poly->Tau[Y] = poly->Tau[Y] + Tau[Y];
                    poly->Tau[Z] = poly->Tau[Z] + Tau[Z];
                }
            }
        }
        /* calibrate the sum of discrete forces into integration by ds */
        poly->Fp[X] = -poly->Fp[X] * poly->area / nodeN;
        poly->Fp[Y] = -poly->Fp[Y] * poly->area / nodeN;
        poly->Fp[Z] = -poly->Fp[Z] * poly->area / nodeN;
        poly->Fv[X] = -poly->Fv[X] * poly->area / nodeN;
        poly->Fv[Y] = -poly->Fv[Y] * poly->area / nodeN;
        poly->Fv[Z] = -poly->Fv[Z] * poly->area / nodeN;
        poly->Tau[X] = -poly->Tau[X] * poly->area / nodeN;
        poly->Tau[Y] = -poly->Tau[Y] * poly->area / nodeN;
        poly->Tau[Z] = -poly->Tau[Z] * poly->area / nodeN;
        return;
    }
}
/* a good practice: end file with a newline */

