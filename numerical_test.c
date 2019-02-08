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
    TCN = 2, /* position index of center node in stencil */
    TTN = 5, /* number of nodes in a stencil */
} TestConst;
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void ComputeSolutionError(Space *space)
{
    FILE *fp = Fopen("solution_error.csv", "w");
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    Real *restrict Us = NULL; /* numerical solution */
    Real *restrict Ue = NULL; /* exact solution */
    int idx = 0; /* linear array index math variable */
    const int meshN = MaxInt(part->m[X], MaxInt(part->m[Y], part->m[Z]));
    Real norm[3] = {0.0}; /* Lp norms */
    int N = 0; /* number of nodes */
    Real err = 0.0; /* solution error */
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
    fprintf(fp, "# mesh, l1 norm, l2 norm, max norm\n");
    fprintf(fp, "%d, %.6g, %.6g, %.6g\n", meshN, norm[1], norm[2], norm[0]);
    fclose(fp);
    return;
}
void ComputeSolutionFunctional(const Time *time, Space *space, const Model *model)
{
    FILE *fp = Fopen("solution_functional.csv", "a");
    if (0 == time->stepC) { /* initialization step */
        fprintf(fp, "# time, kinetic energy, enstrophy \n");
    }
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    Real *restrict U = NULL; /* numerical solution */
    int idx = 0; /* linear array index math variable */
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    Real rho = model->refRho; /* density */
    RealVec V = {0.0}; /* velocity */
    RealVec W = {0.0}; /* vorticity */
    Real Vs[DIMS][TTN] = {{0.0}}; /* velocity stencil */
    Real dV[DIMS][DIMS] = {{0.0}}; /* velocity gradient */
    Real Ek = 0.0; /* kinetic energy */
    Real Ee = 0.0; /* enstrophy */
    int N = 0; /* number of nodes */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                for (int s = 0; s < DIMS; ++s) {
                    for (int n = -TCN; n <= TCN; ++n) {
                        idx = IndexNode(k + n * h[s][Z], j + n * h[s][Y], i + n * h[s][X], part->n[Y], part->n[X]);
                        U = node[idx].U[TO];
                        Vs[X][TCN+n] = U[1] / U[0];
                        Vs[Y][TCN+n] = U[2] / U[0];
                        Vs[Z][TCN+n] = U[3] / U[0];
                    }
                    dV[X][s] = (-Vs[X][TCN+2] + 8.0 * Vs[X][TCN+1] - 8.0 * Vs[X][TCN-1] + Vs[X][TCN-2]) / (12.0 * d[s]);
                    dV[Y][s] = (-Vs[Y][TCN+2] + 8.0 * Vs[Y][TCN+1] - 8.0 * Vs[Y][TCN-1] + Vs[Y][TCN-2]) / (12.0 * d[s]);
                    dV[Z][s] = (-Vs[Z][TCN+2] + 8.0 * Vs[Z][TCN+1] - 8.0 * Vs[Z][TCN-1] + Vs[Z][TCN-2]) / (12.0 * d[s]);
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
    fprintf(fp, "%.6g, %.6g, %.6g\n", time->now, Ek, Ee);
    fclose(fp);
    return;
}
/* a good practice: end file with a newline */

