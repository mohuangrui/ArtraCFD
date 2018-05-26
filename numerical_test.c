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
#include "numerical_test.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include <math.h> /* common mathematical functions */
#include "boundary_treatment.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    CASEN = 2, /* total number of test cases */
    CASE = 1, /* current case index */
    CEN = 2, /* position index of center node in stencil */
    NSTENCIL = 5, /* number of nodes in a stencil */
} CaseConstants;
/****************************************************************************
 * Function Pointers
 ****************************************************************************/
typedef void (*TestCase)(const Real, const Real, const Real, 
        const Model *, const Real [restrict], Real [restrict]);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void VortexPreservation(const Real, const Real, const Real, 
        const Model *, const Real [restrict], Real [restrict]);
static void TaylorGreenVortex(const Real, const Real, const Real, 
        const Model *, const Real [restrict], Real [restrict]);
/****************************************************************************
 * Global Variables Definition with Private Scope
 ****************************************************************************/
static TestCase SetCase[CASEN] = {
    VortexPreservation,
    TaylorGreenVortex};
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int SetField(const int tn, Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    /* extract global initial values */
    const Real Uo[DIMUo] = {
        part->valueIC[0][ENTRYIC-5],
        part->valueIC[0][ENTRYIC-4],
        part->valueIC[0][ENTRYIC-3],
        part->valueIC[0][ENTRYIC-2],
        part->valueIC[0][ENTRYIC-1]};
    Real Ue[DIMUo] = {0.0};
    RealVec p = {0.0};
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                p[X] = PointSpace(i, part->domain[X][MIN], part->d[X], part->ng);
                p[Y] = PointSpace(j, part->domain[Y][MIN], part->d[Y], part->ng);
                p[Z] = PointSpace(k, part->domain[Z][MIN], part->d[Z], part->ng);
                SetCase[CASE](p[X], p[Y], p[Z], model, Uo, Ue);
                ConservativeByPrimitive(model->gamma, Ue, node[idx].U[tn]);
            }
        }
    }
    BoundaryConditionsAndTreatments(tn, space, model);
    return 0;
}
static void VortexPreservation(const Real x, const Real y, const Real z, 
        const Model *model, const Real Uo[restrict], Real Ue[restrict])
{
    /* note: model->gasR should be 1.0 */
    const Real pi = 3.14159265359;
    const Real g = model->gamma;
    const Real G = 5.0; /* vortex strength */
    const Real R = 5.0; /* vortex radius */
    const Real r2 = x * x + y * y + 0.0 * z;
    if (R * R < r2) {
        for (int dim = 0; dim < DIMUo; ++dim) {
            Ue[dim] = Uo[dim];
        }
        return;
    }
    Ue[1] = Uo[1] + (G / (2.0 * pi)) * exp(0.5 * (1.0 - r2)) * (-y);
    Ue[2] = Uo[2] + (G / (2.0 * pi)) * exp(0.5 * (1.0 - r2)) * (x);
    Ue[3] = 0.0;
    Ue[5] = Uo[4] / Uo[0] - (g - 1.0) * G * G / (8.0 * g * pi * pi) * exp(1.0 - r2);
    Ue[0] = pow(Ue[5], 1.0 / (g - 1.0));
    Ue[4] = pow(Ue[5], g / (g - 1.0));
    return;
}
static void TaylorGreenVortex(const Real x, const Real y, const Real z, 
        const Model *model, const Real Uo[restrict], Real Ue[restrict])
{
    /* note: model->refMu should be constant */
    Ue[0] = Uo[0];
    Ue[1] = Uo[1] * sin(x) * cos(y) * cos(z);
    Ue[2] = -Uo[2] * cos(x) * sin(y) * cos(z);
    Ue[3] = Uo[3];
    Ue[4] = Uo[4] + (Uo[0] * Uo[1] * Uo[1] / 16.0) * (cos(2.0*x) + cos(2*y)) * (cos(2*z) + 2);
    Ue[5] = Ue[4] / (Ue[0] * model->gasR);
    return;
}
int ComputeSolutionError(Space *space, const Model *model)
{
    FILE *filePointer = fopen("solution_error.csv", "w");
    if (NULL == filePointer) {
        FatalError("failed to write data...");
    }
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    Real *restrict Us = NULL; /* numerical solution */
    Real *restrict Ue = NULL; /* exact solution */
    int idx = 0; /* linear array index math variable */
    const int meshN = MaxInt(part->m[X], MaxInt(part->m[Y], part->m[Z]));
    Real norm[3] = {0.0}; /* Lp norms */
    int N = 0; /* number of nodes */
    Real err = 0.0; /* solution error */
    SetField(TN, space, model); /* compute exact solution field */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                Us = node[idx].U[TO];
                Ue = node[idx].U[TN];
                err = fabs(Us[0] - Ue[0]);
                norm[0] = MaxReal(norm[0], err);
                norm[1] = norm[1] + err;
                norm[2] = norm[2] + err * err;
                ++N;
            }
        }
    }
    norm[1] = norm[1] / N;
    norm[2] = sqrt(norm[2] / N);
    fprintf(filePointer, "# mesh, l1 norm, l2 norm, max norm\n"); 
    fprintf(filePointer, "%d, %.6g, %.6g, %.6g\n", meshN, norm[1], norm[2], norm[0]); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
int ComputeSolutionFunctional(const Time *time, Space *space, const Model *model)
{
    FILE *filePointer = fopen("solution_functional.csv", "a");
    if (NULL == filePointer) {
        FatalError("failed to write data...");
    }
    if (0 == time->stepC) { /* this is the initialization step */
        fprintf(filePointer, "# time, kinetic energy, enstrophy \n"); 
    }
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    Real *restrict U = NULL; /* numerical solution */
    int idx = 0; /* linear array index math variable */
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    Real rho = model->refRho; /* density */
    RealVec V = {0.0}; /* velocity */
    RealVec W = {0.0}; /* vorticity */
    Real Vs[DIMS][NSTENCIL] = {{0.0}}; /* velocity stencil */
    Real dV[DIMS][DIMS] = {{0.0}}; /* velocity gradient */
    Real Ek = 0.0; /* kinetic energy */
    Real Ee = 0.0; /* enstrophy */
    int N = 0; /* number of nodes */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                for (int s = 0; s < DIMS; ++s) {
                    for (int n = -CEN; n <= CEN; ++n) {
                        idx = IndexNode(k + n * h[s][Z], j + n * h[s][Y], i + n * h[s][X], part->n[Y], part->n[X]);
                        U = node[idx].U[TO];
                        Vs[X][CEN+n] = U[1] / U[0];
                        Vs[Y][CEN+n] = U[2] / U[0];
                        Vs[Z][CEN+n] = U[3] / U[0];
                    }
                    dV[X][s] = (-Vs[X][CEN+2] + 8.0 * Vs[X][CEN+1] - 8.0 * Vs[X][CEN-1] + Vs[X][CEN-2]) / (12.0 * d[s]);
                    dV[Y][s] = (-Vs[Y][CEN+2] + 8.0 * Vs[Y][CEN+1] - 8.0 * Vs[Y][CEN-1] + Vs[Y][CEN-2]) / (12.0 * d[s]);
                    dV[Z][s] = (-Vs[Z][CEN+2] + 8.0 * Vs[Z][CEN+1] - 8.0 * Vs[Z][CEN-1] + Vs[Z][CEN-2]) / (12.0 * d[s]);
                }
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                U = node[idx].U[TO];
                rho = U[0];
                V[X] = U[1] / U[0];
                V[Y] = U[2] / U[0];
                V[Z] = U[3] / U[0];
                W[X] = dV[Z][Y] - dV[Y][Z];
                W[Y] = dV[X][Z] - dV[Z][X];
                W[Z] = dV[Y][X] - dV[X][Y];
                Ek = Ek + rho * Dot(V,V);
                Ee = Ee + rho * Dot(W,W);
                ++N;
            }
        }
    }
    Ek = 0.5 * Ek / N;
    Ee = 0.5 * Ee / N;
    fprintf(filePointer, "%.6g, %.6g, %.6g\n", time->now, Ek, Ee); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

