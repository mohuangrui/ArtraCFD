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
#include "convective_flux.h"
#include "weno.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    FDN = 5, /* width of the direct stencil */
    FTN = 6, /* width of the entire stencil */
} FhatConst;
/****************************************************************************
 * Function Pointers
 ****************************************************************************/
typedef void (*FhatReconstructor)(Real [restrict][DIMU], Real [restrict]);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void CharacteristicVariable(const int, const int, const int, const int,
        const int, const int, const int, const int [restrict], const Node *const,
        Real [restrict][DIMU], Real [restrict][DIMU]);
static void CharacteristicFlux(const Real [restrict], Real [restrict][DIMU],
        const int, const int, const int,  Real [restrict][DIMU]);
static void InverseProjection(Real [restrict][DIMU], const Real [restrict],
        const Real [restrict], Real [restrict]);
/****************************************************************************
 * Global Variables Definition with Private Scope
 ****************************************************************************/
static FhatReconstructor ReconstructFhat[2] = {
    WENO3,
    WENO5};
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void ComputeFhat(const int tn, const int s, const int k, const int j, const int i,
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
    Real W[FTN][DIMU];
    CharacteristicVariable(tn, s, k, j, i, model->sL, model->sR, partn, node, L, W);
    /* construct local characteristic fluxes */
    Real HP[FDN][DIMU]; /* forward characteristic flux stencil */
    Real HN[FDN][DIMU]; /* backward characteristic flux stencil */
    CharacteristicFlux(LambdaP, W, 0, +1, model->sR - model->sL, HP);
    CharacteristicFlux(LambdaN, W, model->sR - model->sL, -1, model->sR - model->sL, HN);
    /* WENO reconstruction */
    Real HhatP[DIMU]; /* forward numerical flux of characteristic fields */
    Real HhatN[DIMU]; /* backward numerical flux of characteristic fields */
    ReconstructFhat[model->sScheme](HP, HhatP);
    ReconstructFhat[model->sScheme](HN, HhatN);
    /* inverse projection */
    InverseProjection(R, HhatP, HhatN, Fhat);
    return;
}
static void CharacteristicVariable(const int tn, const int s, const int k, const int j,
        const int i, const int sL, const int sR, const int partn[restrict],
        const Node *const node, Real L[restrict][DIMU], Real W[restrict][DIMU])
{
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    int idx = 0; /* linear array index math variable */
    const Real *restrict U = NULL;
    for (int n = sL, m = 0; n <= sR; ++n, ++m) {
        idx = IndexNode(k + n * h[s][Z], j + n * h[s][Y], i + n * h[s][X], partn[Y], partn[X]);
        U = node[idx].U[tn];
        for (int r = 0; r < DIMU; ++r) {
            W[m][r] = 0.0;
            for (int c = 0; c < DIMU; ++c) {
                W[m][r] = W[m][r] + L[r][c] * U[c];
            }
        }
    }
    return;
}
static void CharacteristicFlux(const Real Lambda[restrict], Real W[restrict][DIMU],
        const int start, const int wind, const int tot, Real H[restrict][DIMU])
{
    for (int n = start, m = 0; m < tot; n = n + wind, ++m) {
        for (int r = 0; r < DIMU; ++r) {
            H[m][r] = Lambda[r] * W[n][r];
        }
    }
    return;
}
static void InverseProjection(Real R[restrict][DIMU], const Real HhatP[restrict],
        const Real HhatN[restrict], Real Fhat[restrict])
{
    for (int r = 0; r < DIMU; ++r) {
        Fhat[r] = 0.0;
        for (int c = 0; c < DIMU; ++c) {
            Fhat[r] = Fhat[r] + R[r][c] * (HhatP[c] + HhatN[c]);
        }
    }
    return;
}
/* a good practice: end file with a newline */

