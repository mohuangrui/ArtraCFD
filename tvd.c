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
#include "tvd.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int CalculateReconstructedFlux(
        Real [], const Real [], const Real [], Real [][DIMU], const Real []);
static int FluxDecompositionCoefficientPhi(
        const int, Real [], const Real, const int, const int, const int, 
        const Real *, const Space *, const Model *);
static int FunctionG(
        const int, Real [], const Real, const int, const int, const int, 
        const Real *, const Space *, const Model *);
static int CalculateGamma(Real [], const Real [], const Real [], const Real []);
static int CalculateSigma(Real [], const Real [], const Real, const Real);
static Real NumericalDissipationDelta(const Real [], const Real);
static Real Q(const Real, const Real);
static Real minmod(const Real, const Real);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int TVD(const int s, Real Fhat[], const Real r, const int k, const int j,
        const int i, const Real *U, const Space *space, const Model *model)
{
    Real F[DIMU] = {0.0}; /* flux at current node */
    Real Fh[DIMU] = {0.0}; /* flux at neighbour */
    Real R[DIMU][DIMU] = {{0.0}}; /* vector space {Rn} */
    Real Phi[DIMU] = {0.0}; /* flux projection or decomposition coefficients on {Rn} */
    Real Uo[DIMUo] = {0.0}; /* averaged rho, u, v, w, hT, c */
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxh = IndexMath(k + h[s][Z], j + h[s][Y], i + h[s][X], space) * DIMU;
    ConvectiveFlux(s, F, idx, U, model->gamma);
    ConvectiveFlux(s, Fh, idxh, U, model->gamma);
    SymmetricAverage(Uo, idx, idxh, U, model->gamma, model->averager);
    EigenvectorSpaceR(s, R, Uo);
    FluxDecompositionCoefficientPhi(s, Phi, r, k, j, i, U, space, model);
    CalculateReconstructedFlux(Fhat, F, Fh, R, Phi);
    return 0;
}
static int CalculateReconstructedFlux(
        Real Fhat[], const Real F[], const Real Fh[], Real R[][DIMU], const Real Phi[])
{
    Real RPhi[DIMU] = {0.0}; /* R x Phi */
    for (int row = 0; row < DIMU; ++row) {
        RPhi[row] = 0;
        for (int dummy = 0; dummy < DIMU; ++dummy) {
            RPhi[row] = RPhi[row] + R[row][dummy] * Phi[dummy];
        }
    }
    for (int row = 0; row < DIMU; ++row) {
        Fhat[row] = 0.5 * (F[row] + Fh[row] + RPhi[row]);
    }
    return 0;
}
static int FluxDecompositionCoefficientPhi(const int s, Real Phi[], const Real r,
        const int k, const int j, const int i, const Real *U, const Space *space, const Model *model)
{
    Real g[DIMU] = {0.0}; /* TVD function g at current node */
    Real gh[DIMU] = {0.0}; /* TVD function g at neighbour */
    Real gamma[DIMU] = {0.0}; /* TVD function gamma */
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real L[DIMU][DIMU] = {{0.0}}; /* vector space {Ln} */
    Real alpha[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real Uo[DIMUo] = {0.0}; /* averaged rho, u, v, w, hT, c */
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxh = IndexMath(k + h[s][Z], j + h[s][Y], i + h[s][X], space) * DIMU;
    SymmetricAverage(Uo, idx, idxh, U, model->gamma, model->averager);
    EigenvalueLambda(s, lambda, Uo);
    EigenvectorSpaceL(s, L, Uo, model->gamma);
    DecompositionCoefficientAlpha(alpha, L, idx, idxh, U);
    FunctionG(s, g, r, k, j, i, U, space, model);
    FunctionG(s, gh, r, k + h[s][Z], j + h[s][Y], i + h[s][X], U, space, model);
    const Real delta = NumericalDissipationDelta(Uo, model->delta); /* numerical dissipation */
    CalculateGamma(gamma, g, gh, alpha);
    for (int row = 0; row < DIMU; ++row) {
        Phi[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta) * alpha[row];
    }
    return 0;
}
static int FunctionG(const int s, Real g[], const Real r, const int k, const int j, const int i,
        const Real *U, const Space *space, const Model *model)
{
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real lambdah[DIMU] = {0.0}; /* eigenvalues at neighbour */
    Real L[DIMU][DIMU] = {{0.0}}; /* vector space {Ln} */
    Real Lh[DIMU][DIMU] = {{0.0}}; /* vector space {Ln} */
    Real alpha[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real alphah[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real sigma[DIMU] = {0.0}; /* TVD function sigma */
    Real sigmah[DIMU] = {0.0}; /* TVD function sigma at neighbour */
    Real Uo[DIMUo] = {0.0}; /* averaged rho, u, v, w, hT, c */
    Real Uoh[DIMUo] = {0.0}; /* averaged rho, u, v, w, hT, c */
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    const int idxl = IndexMath(k - h[s][Z], j - h[s][Y], i - h[s][X], space) * DIMU;
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxr = IndexMath(k + h[s][Z], j + h[s][Y], i + h[s][X], space) * DIMU;
    SymmetricAverage(Uo, idx, idxr, U, model->gamma, model->averager);
    SymmetricAverage(Uoh, idxl, idx, U, model->gamma, model->averager);
    EigenvalueLambda(s, lambda, Uo);
    EigenvalueLambda(s, lambdah, Uoh);
    EigenvectorSpaceL(s, L, Uo, model->gamma);
    EigenvectorSpaceL(s, Lh, Uoh, model->gamma);
    DecompositionCoefficientAlpha(alpha, L, idx, idxr, U);
    DecompositionCoefficientAlpha(alphah, Lh, idxl, idx, U);
    const Real delta = NumericalDissipationDelta(Uo, model->delta); /* numerical dissipation */
    const Real deltah = NumericalDissipationDelta(Uoh, model->delta); /* numerical dissipation */
    CalculateSigma(sigma, lambda, r, delta);
    CalculateSigma(sigmah, lambdah, r, deltah);
    for (int row = 0; row < DIMU; ++row) {
        g[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
    }
    return 0;
}
static int CalculateGamma(Real gamma[], const Real g[], const Real gh[], const Real alpha[])
{
    for (int row = 0; row < DIMU; ++row) {
        if (0 != alpha[row]) {
            gamma[row] = (gh[row] - g[row]) / alpha[row];
        }
        else {
            gamma[row] = 0;
        }
    }
    return 0;
}
static int CalculateSigma(Real sigma[], const Real lambda[], const Real r, const Real delta)
{
    for (int row = 0; row < DIMU; ++row) {
        sigma[row] = 0.5 * (Q(lambda[row], delta) - r * lambda[row] * lambda[row]);
    }
    return 0;
}
static Real NumericalDissipationDelta(const Real Uo[], const Real delta0)
{
    return delta0 * (fabs(Uo[1]) + fabs(Uo[2]) + fabs(Uo[3]) + Uo[5]);
}
static Real Q(const Real z, const Real delta)
{
    if (fabs(z) >= delta) {
        return fabs(z);
    }
    return (0.5 * (z * z / delta + delta));
}
static Real minmod(const Real x, const Real y)
{
    return (Sign(x) * MaxReal(0, MinReal(fabs(x), y * Sign(x))));
}
/* a good practice: end file with a newline */

