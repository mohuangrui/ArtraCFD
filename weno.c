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
        const int [restrict], const Node *const , Real [restrict][DIMU], const Real [restrict],
        const int, const int, Real [restrict][NSTENCIL]);
static void WENO5(Real [restrict][NSTENCIL], Real [restrict]);
static void InverseProjection(Real [restrict][DIMU], const Real [restrict], 
        const Real [restrict], Real [restrict]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void WENO(const int tn, const int s, const int k, const int j, const int i, 
        const int partn[restrict], const Node *const node, const Model *model, Real Fhat[restrict])
{
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    const int idxL = IndexNode(k, j, i, partn[Y], partn[X]);
    const int idxR = IndexNode(k + h[s][Z], j + h[s][Y], i + h[s][X], partn[Y], partn[X]);
    /* evaluate interface values by averaging */
    Real Uo[DIMUo] = {0.0}; /* store averaged primitives */
    SymmetricAverage(model->averager, model->gamma, node[idxL].U[tn], node[idxR].U[tn], Uo);
    /* decompose Jocabian matrix */
    Real Lambda[DIMU] = {0.0}; /* eigenvalues */
    Real L[DIMU][DIMU] = {{0.0}}; /* vector space {Ln} */
    Real R[DIMU][DIMU] = {{0.0}}; /* vector space {Rn} */
    Eigenvalue(s, Uo, Lambda);
    EigenvectorL(s, model->gamma, Uo, L);
    EigenvectorR(s, Uo, R);
    /* flux vector splitting */
    Real LambdaP[DIMU] = {0.0}; /* eigenvalues */
    Real LambdaN[DIMU] = {0.0}; /* eigenvalues */
    EigenvalueSplitting(model->splitter, Lambda, LambdaP, LambdaN);
    /* construct local characteristic fluxes */
    Real HP[DIMU][NSTENCIL] = {{0.0}}; /* forward characteristic flux stencil */
    Real HN[DIMU][NSTENCIL] = {{0.0}}; /* backward characteristic flux stencil */
    CharacteristicProjection(tn, s, k, j, i, partn, node, L, LambdaP, -2, +1, HP);
    CharacteristicProjection(tn, s, k, j, i, partn, node, L, LambdaN, +3, -1, HN);
    /* WENO reconstruction */
    Real HhatP[DIMU] = {0.0}; /* forward numerical flux of characteristic fields */
    Real HhatN[DIMU] = {0.0}; /* backward numerical flux of characteristic fields */
    WENO5(HP, HhatP);
    WENO5(HN, HhatN);
    /* inverse projection */
    InverseProjection(R, HhatP, HhatN, Fhat);
    return;
}
static void CharacteristicProjection(const int tn, const int s, const int k, const int j, const int i, 
        const int partn[restrict], const Node *const node, Real L[restrict][DIMU], const Real Lambda[restrict],
        const int startN, const int wind, Real H[restrict][NSTENCIL])
{
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    int idx = 0; /* linear array index math variable */
    const Real *restrict U = NULL;
    for (int n = startN, count = 0; count < NSTENCIL; n = n + wind, ++count) {
        idx = IndexNode(k + n * h[s][Z], j + n * h[s][Y], i + n * h[s][X], partn[Y], partn[X]);
        U = node[idx].U[tn];
        for (int row = 0; row < DIMU; ++row) {
            H[row][count] = 0.0;
            for (int dummy = 0; dummy < DIMU; ++dummy) {
                H[row][count] = H[row][count] + L[row][dummy] * U[dummy];
            }
            H[row][count] = H[row][count] * Lambda[row];
        }
    }
    return;
}
static void WENO5(Real F[restrict][NSTENCIL], Real Fhat[restrict])
{
    Real omega[R] = {0.0}; /* weights */
    Real q[R] = {0.0}; /* q vectors */
    Real alpha[R] = {0.0};
    Real IS[R] = {0.0};
    const Real C[R] = {0.1, 0.6, 0.3};
    const Real epsilon = 1.0e-6;
    for (int row = 0; row < DIMU; ++row) {
        IS[0] = (13.0 / 12.0) * Square(F[row][CEN-2] - 2.0 * F[row][CEN-1] + F[row][CEN]) + 
            (1.0 / 4.0) * Square(F[row][CEN-2] - 4.0 * F[row][CEN-1] + 3.0 * F[row][CEN]);
        IS[1] = (13.0 / 12.0) * Square(F[row][CEN-1] - 2.0 * F[row][CEN] + F[row][CEN+1]) +
            (1.0 / 4.0) * Square(F[row][CEN-1] - F[row][CEN+1]);
        IS[2] = (13.0 / 12.0) * Square(F[row][CEN] - 2.0 * F[row][CEN+1] + F[row][CEN+2]) +
            (1.0 / 4.0) * Square(3.0 * F[row][CEN] - 4.0 * F[row][CEN+1] + F[row][CEN+2]);
        alpha[0] = C[0] / Square(epsilon + IS[0]);
        alpha[1] = C[1] / Square(epsilon + IS[1]);
        alpha[2] = C[2] / Square(epsilon + IS[2]);
        omega[0] = alpha[0] / (alpha[0] + alpha[1] + alpha[2]);
        omega[1] = alpha[1] / (alpha[0] + alpha[1] + alpha[2]);
        omega[2] = alpha[2] / (alpha[0] + alpha[1] + alpha[2]);
        q[0] = (1.0 / 6.0) * (2.0 * F[row][CEN-2] - 7.0 * F[row][CEN-1] + 11.0 * F[row][CEN]);
        q[1] = (1.0 / 6.0) * (-F[row][CEN-1] + 5.0 * F[row][CEN] + 2.0 * F[row][CEN+1]);
        q[2] = (1.0 / 6.0) * (2.0 * F[row][CEN] + 5.0 * F[row][CEN+1] - F[row][CEN+2]);
        Fhat[row] = omega[0] * q[0] + omega[1] * q[1] + omega[2] * q[2];
    }
    return;
}
static void InverseProjection(Real R[restrict][DIMU], const Real HhatP[restrict], 
        const Real HhatN[restrict], Real Fhat[restrict])
{
    for (int row = 0; row < DIMU; ++row) {
        Fhat[row] = 0.0;
        for (int dummy = 0; dummy < DIMU; ++dummy) {
            Fhat[row] = Fhat[row] + R[row][dummy] * (HhatP[dummy] + HhatN[dummy]);
        }
    }
    return;
}
/* a good practice: end file with a newline */

