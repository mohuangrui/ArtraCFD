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
    R = 2, /* WENO r */
    N = 1, /* N = (r - 1) */
    CEN = 1, /* position index of center node in stencil */
    NSTENCIL = 3, /* number of nodes in a stencil = (2r - 1) */
    TNSTENCIL = 4, /* total number of involved nodes = (2r) */
} WENOConstants;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void CharacteristicVariable(const int, const int, const int, const int, const int, 
        const int [restrict], const Node *const, Real [restrict][DIMU], Real [restrict][DIMU]);
static void CharacteristicFlux(const Real [restrict], Real [restrict][DIMU],
        const int, const int, Real [restrict][DIMU]);
static void WENOConstruction(Real [restrict][DIMU], Real [restrict]);
static void InverseProjection(Real [restrict][DIMU], const Real [restrict], 
        const Real [restrict], Real [restrict]);
static Real Square(const Real);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Jiang, G.S. and Shu, C.W., 1996. Efficient Implementation of Weighted
 * ENO Schemes. Journal of Computational Physics, 126(1), pp.202-228.
 */
void WENO3(const int tn, const int s, const int k, const int j, const int i, 
        const int partn[restrict], const Node *const node, const Model *model, Real Fhat[restrict])
{
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    const int idxL = IndexNode(k, j, i, partn[Y], partn[X]);
    const int idxR = IndexNode(k + h[s][Z], j + h[s][Y], i + h[s][X], partn[Y], partn[X]);
    /* evaluate interface values by averaging */
    Real Uo[DIMUo]; /* store averaged primitives */
    SymmetricAverage(model->jacobMean, model->gamma, node[idxL].U[tn], node[idxR].U[tn], Uo);
    /* decompose Jacobian matrix */
    Real Lambda[DIMU]; /* eigenvalues */
    Real L[DIMU][DIMU]; /* vector space {Ln} */
    Real R[DIMU][DIMU]; /* vector space {Rn} */
    Eigenvalue(s, Uo, Lambda);
    EigenvectorL(s, model->gamma, Uo, L);
    EigenvectorR(s, Uo, R);
    /* flux vector splitting */
    Real LambdaP[DIMU]; /* eigenvalues */
    Real LambdaN[DIMU]; /* eigenvalues */
    EigenvalueSplitting(model->fluxSplit, Lambda, LambdaP, LambdaN);
    /* construct local characteristic variables for all potential stencils */
    Real W[TNSTENCIL][DIMU];
    CharacteristicVariable(tn, s, k, j, i, partn, node, L, W);
    /* construct local characteristic fluxes */
    Real HP[NSTENCIL][DIMU]; /* forward characteristic flux stencil */
    Real HN[NSTENCIL][DIMU]; /* backward characteristic flux stencil */
    CharacteristicFlux(LambdaP, W, 0, +1, HP);
    CharacteristicFlux(LambdaN, W, NSTENCIL, -1, HN);
    /* WENO reconstruction */
    Real HhatP[DIMU]; /* forward numerical flux of characteristic fields */
    Real HhatN[DIMU]; /* backward numerical flux of characteristic fields */
    WENOConstruction(HP, HhatP);
    WENOConstruction(HN, HhatN);
    /* inverse projection */
    InverseProjection(R, HhatP, HhatN, Fhat);
    return;
}
static void CharacteristicVariable(const int tn, const int s, const int k, const int j, const int i, 
        const int partn[restrict], const Node *const node, Real L[restrict][DIMU], Real W[restrict][DIMU])
{
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    int idx = 0; /* linear array index math variable */
    const Real *restrict U = NULL;
    for (int n = -N, count = 0; count < TNSTENCIL; ++n, ++count) {
        idx = IndexNode(k + n * h[s][Z], j + n * h[s][Y], i + n * h[s][X], partn[Y], partn[X]);
        U = node[idx].U[tn];
        for (int row = 0; row < DIMU; ++row) {
            W[count][row] = 0.0;
            for (int dummy = 0; dummy < DIMU; ++dummy) {
                W[count][row] = W[count][row] + L[row][dummy] * U[dummy];
            }
        }
    }
    return;
}
static void CharacteristicFlux(const Real Lambda[restrict], Real W[restrict][DIMU],
        const int startN, const int wind, Real H[restrict][DIMU])
{
    for (int n = startN, count = 0; count < NSTENCIL; n = n + wind, ++count) {
        for (int row = 0; row < DIMU; ++row) {
            H[count][row] = Lambda[row] * W[n][row];
        }
    }
    return;
}
static void WENOConstruction(Real F[restrict][DIMU], Real Fhat[restrict])
{
    Real omega[R]; /* weights */
    Real q[R]; /* q vectors */
    Real alpha[R];
    Real IS[R];
    const Real C[R] = {1.0 / 3.0, 2.0 / 3.0};
    const Real epsilon = 1.0e-6;
    for (int row = 0; row < DIMU; ++row) {
        IS[0] = Square(F[CEN][row] - F[CEN-1][row]);
        IS[1] = Square(F[CEN+1][row] - F[CEN][row]);
        alpha[0] = C[0] / Square(epsilon + IS[0]);
        alpha[1] = C[1] / Square(epsilon + IS[1]);
        omega[0] = alpha[0] / (alpha[0] + alpha[1]);
        omega[1] = alpha[1] / (alpha[0] + alpha[1]);
        q[0] = (1.0 / 2.0) * (-F[CEN-1][row] + 3.0 * F[CEN][row]);
        q[1] = (1.0 / 2.0) * (F[CEN][row] + F[CEN+1][row]);
        Fhat[row] = omega[0] * q[0] + omega[1] * q[1];
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
static Real Square(const Real x)
{
    return x * x;
}
/* a good practice: end file with a newline */

