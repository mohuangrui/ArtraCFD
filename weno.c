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
#include "weno.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    R = 3, /* WENO r */
    N = 2, /* N = (r - 1) */
    CEN = 2, /* position index of center node in stencil */
    NSTENCIL = 5, /* number of nodes in a stencil = (2r - 1) */
} WENOConstants;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void CharacteristicProjection(const int, const int, const int, const int, const int, 
        const int [restrict], const Node *, Real [restrict][DIMU], const Real [restrict],
        const int , const int , Real [restrict][DIMU]);
static void WENO5(Real [restrict][DIMU], Real [restrict]);
static void InverseProjection(Real [restrict][DIMU], const Real [restrict], 
        const Real [restrict], Real [restrict]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void WENO(const int s, const int tn, const int k, const int j, const int i, 
        const int partn[restrict], const Node *node, const Model *model, Real Fhat[restrict])
{
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    const int idxL = IndexNode(k, j, i, partn[Y], partn[X]);
    const int idxR = IndexNode(k + h[s][Z], j + h[s][Y], i + h[s][X], partn[Y], partn[X]);
    /* evaluate interface values by averaging */
    Real Uo[DIMUo] = {0.0}; /* store averaged primitives */
    SymmetricAverage(model->averager, model->gamma, node[idxL].U[tn], node[idxR].U[tn], Uo);
    /* decompose Jocabian matrix */
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real L[DIMU][DIMU] = {{0.0}}; /* vector space {Ln} */
    Real R[DIMU][DIMU] = {{0.0}}; /* vector space {Rn} */
    Eigenvalue(s, Uo, lambda);
    EigenvectorL(s, model->gamma, Uo, L);
    EigenvectorR(s, Uo, R);
    /* flux vector splitting */
    Real lambdaP[DIMU] = {0.0}; /* eigenvalues */
    Real lambdaN[DIMU] = {0.0}; /* eigenvalues */
    EigenvalueSplitting(model->splitter, lambda, lambdaP, lambdaN);
    /* construct local characteristic fluxes */
    Real HP[NSTENCIL][DIMU] = {{0.0}}; /* forward characteristic flux stencil */
    Real HN[NSTENCIL][DIMU] = {{0.0}}; /* backward characteristic flux stencil */
    CharacteristicProjection(s, tn, k, j, i, partn, node, L, lambdaP, -2, +1, HP);
    CharacteristicProjection(s, tn, k, j, i, partn, node, L, lambdaN, +3, -1, HN);
    Real HhatP[DIMU] = {0.0}; /* forward numerical flux of characteristic fields */
    Real HhatN[DIMU] = {0.0}; /* backward numerical flux of characteristic fields */
    /* WENO reconstruction */
    WENO5(HP, HhatP);
    WENO5(HN, HhatN);
    /* inverse projection */
    InverseProjection(R, HhatP, HhatN, Fhat);
    return;
}
static void CharacteristicProjection(const int s, const int tn, const int k, const int j, const int i, 
        const int partn[restrict], const Node *node, Real L[restrict][DIMU], const Real lambda[restrict],
        const int startN, const int wind, Real H[restrict][DIMU])
{
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    int idx = 0; /* linear array index math variable */
    const Real *restrict U = NULL;
    for (int n = startN, count = 0; count < NSTENCIL; n = n + wind, ++count) {
        idx = IndexNode(k + n * h[s][Z], j + n * h[s][Y], i + n * h[s][X], partn[Y], partn[X]);
        U = node[idx].U[tn];
        for (int row = 0; row < DIMU; ++row) {
            H[count][row] = 0;
            for (int dummy = 0; dummy < DIMU; ++dummy) {
                H[count][row] = H[count][row] + L[row][dummy] * U[dummy];
            }
            H[count][row] = H[count][row] * lambda[row];
        }
    }
    return;
}
static void WENO5(Real F[restrict][DIMU], Real Fhat[restrict])
{
    Real omega[R] = {0.0}; /* weights */
    Real q[R] = {0.0}; /* q vectors */
    Real alpha[R] = {0.0};
    Real IS[R] = {0.0};
    const Real C[R] = {0.1, 0.6, 0.3};
    const Real epsilon = 1.0e-6;
    for (int row = 0; row < DIMU; ++row) {
        IS[0] = (13.0 / 12.0) * pow((F[CEN-2][row] - 2.0 * F[CEN-1][row] + F[CEN][row]), 2.0) + 
            (1.0 / 4.0) * pow((F[CEN-2][row] - 4.0 * F[CEN-1][row] + 3.0 * F[CEN][row]), 2.0);
        IS[1] = (13.0 / 12.0) * pow((F[CEN-1][row] - 2.0 * F[CEN][row] + F[CEN+1][row]), 2.0) +
            (1.0 / 4.0) * pow((F[CEN-1][row] - F[CEN+1][row]), 2.0);
        IS[2] = (13.0 / 12.0) * pow((F[CEN][row] - 2.0 * F[CEN+1][row] + F[CEN+2][row]), 2.0) +
            (1.0 / 4.0) * pow((3.0 * F[CEN][row] - 4.0 * F[CEN+1][row] + F[CEN+2][row]), 2.0);
        alpha[0] = C[0] / pow((epsilon + IS[0]), 2.0);
        alpha[1] = C[1] / pow((epsilon + IS[1]), 2.0);
        alpha[2] = C[2] / pow((epsilon + IS[2]), 2.0);
        omega[0] = alpha[0] / (alpha[0] + alpha[1] + alpha[2]);
        omega[1] = alpha[1] / (alpha[0] + alpha[1] + alpha[2]);
        omega[2] = alpha[2] / (alpha[0] + alpha[1] + alpha[2]);
        q[0] = (2.0 * F[CEN-2][row] - 7.0 * F[CEN-1][row] + 11.0 * F[CEN][row]) / 6.0;
        q[1] = (-F[CEN-1][row] + 5.0 * F[CEN][row] + 2.0 * F[CEN+1][row]) / 6.0;
        q[2] = (2.0 * F[CEN][row] + 5.0 * F[CEN+1][row] - F[CEN+2][row]) / 6.0;
        Fhat[row] = omega[0] * q[0] + omega[1] * q[1] + omega[2] * q[2];
    }
    return;
}
static void InverseProjection(Real R[restrict][DIMU], const Real HhatP[restrict], 
        const Real HhatN[restrict], Real Fhat[restrict])
{
    for (int row = 0; row < DIMU; ++row) {
        Fhat[row] = 0;
        for (int dummy = 0; dummy < DIMU; ++dummy) {
            Fhat[row] = Fhat[row] + R[row][dummy] * (HhatP[dummy] + HhatN[dummy]);
        }
    }
    return;
}
/* a good practice: end file with a newline */

