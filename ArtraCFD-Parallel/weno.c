/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) 2014-2018 Huangrui Mo <huangrui.mo@gmail.com>              *
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
static int CharacteristicProjection(const int, Real [][DIMU], const Real [], Real [][DIMU],
        const int, const int, const int, const int, const int, const Real *, const Space *);
static int NumericalFlux(Real [], Real [][DIMU]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WENO(const int s, Real Fhat[], const Real r, const int k, const int j,
        const int i, const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* Averaged rho, u, v, w, hT, c */
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real L[DIMU][DIMU] = {{0.0}}; /* vector space {Ln} */
    Real R[DIMU][DIMU] = {{0.0}}; /* vector space {Rn} */
    Real HhatPlus[DIMU] = {0.0}; /* forward numerical flux of characteristic fields */
    Real HhatMinus[DIMU] = {0.0}; /* backward numerical flux of characteristic fields */
    Real HPlus[NSTENCIL][DIMU] = {{0.0}}; /* forward characteristic fields stencil */
    Real HMinus[NSTENCIL][DIMU] = {{0.0}}; /* backward characteristic fields stencil */
    Real lambdaPlus[DIMU] = {0.0}; /* eigenvalues */
    Real lambdaMinus[DIMU] = {0.0}; /* eigenvalues */
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxh = IndexMath(k + h[s][Z], j + h[s][Y], i + h[s][X], space) * DIMU;
    SymmetricAverage(Uo, idx, idxh, U, model->gamma, model->averager);
    EigenvalueLambda(s, lambda, Uo);
    EigenvectorSpaceL(s, L, Uo, model->gamma);
    EigenvectorSpaceR(s, R, Uo);
    FluxSplitting(lambdaPlus, lambdaMinus, lambda, model->splitter);
    CharacteristicProjection(s, HPlus, lambdaPlus, L, -2, +1, k, j, i, U, space);
    CharacteristicProjection(s, HMinus, lambdaMinus, L, +3, -1, k, j, i, U, space);
    NumericalFlux(HhatPlus, HPlus);
    NumericalFlux(HhatMinus, HMinus);
    for (int row = 0; row < DIMU; ++row) {
        Fhat[row] = 0;
        for (int dummy = 0; dummy < DIMU; ++dummy) {
            Fhat[row] = Fhat[row] + R[row][dummy] * (HhatPlus[dummy] + HhatMinus[dummy]);
        }
    }
    return (!r); /* to avoid unsed parameter warning */
}
static int CharacteristicProjection(const int s, Real H[][DIMU], const Real lambda[], Real L[][DIMU],
        const int startN, const int wind, const int k, const int j, const int i, const Real *U, const Space *space)
{
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    int idx = 0;
    for (int n = startN, count = 0; count < NSTENCIL; n = n + wind, ++count) {
        idx = IndexMath(k + n * h[s][Z], j + n * h[s][Y], i + n * h[s][X], space) * DIMU;
        for (int row = 0; row < DIMU; ++row) {
            H[count][row] = 0;
            for (int dummy = 0; dummy < DIMU; ++dummy) {
                H[count][row] = H[count][row] + L[row][dummy] * U[idx + dummy];
            }
            H[count][row] = H[count][row] * lambda[row];
        }
    }
    return 0;
}
static int NumericalFlux(Real Fhat[], Real F[][DIMU])
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
    return 0;
}
/* a good practice: end file with a newline */

