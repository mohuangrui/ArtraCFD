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
#include "cfd_commons.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int CalculateReconstructedFlux(
        Real Fhat[], const Real F[], const Real Fh[], Real R[][DIMU], const Real Phi[]);
static int ComputeFluxDecompositionCoefficientPhiZ(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt);
static int ComputeFluxDecompositionCoefficientPhiY(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt);
static int ComputeFluxDecompositionCoefficientPhiX(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt);
static int ComputeFunctionGZ(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt);
static int ComputeFunctionGY(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt);
static int ComputeFunctionGX(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt);
static int CalculateGamma(
        Real gamma[], const Real g[], const Real gh[], const Real alpha[]);
static int CalculateSigma(
        Real sigma[], const Real lambda[], const Real delta[], const Real r);
static int ComputeNumericalDissipationDeltaZ(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Model *model);
static int ComputeNumericalDissipationDeltaY(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Model *model);
static int ComputeNumericalDissipationDeltaX(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Model *model);
static Real Q(const Real z, const Real delta);
static Real minmod(const Real x, const Real y);
static int sign(const Real x);
static Real min(const Real x, const Real y);
static Real max(const Real x, const Real y);
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaY(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaX(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int CalculateAlpha(
        Real alpha[], Real L[][DIMU], const Real deltaU[]);
static int ComputeEigenvaluesAndEigenvectorSpaceLZ(
        Real lambda[], Real L[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeEigenvaluesAndEigenvectorSpaceLY(
        Real lambda[], Real L[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeEigenvaluesAndEigenvectorSpaceLX(
        Real lambda[], Real L[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeEigenvectorSpaceRZ(
        Real R[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeEigenvectorSpaceRY(
        Real R[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeEigenvectorSpaceRX(
        Real R[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeRoeAverageZ(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeRoeAverageY(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeRoeAverageX(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int CalculateRoeAverageUo(
        Real Uo[], const int idx, const int idxh, 
        const Real *U, const Model *model);
static int ComputeNonViscousFluxZ(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeNonViscousFluxY(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
static int ComputeNonViscousFluxX(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int TVDFluxZ(Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    Real F[DIMU] = {0.0}; /* flux at current node */
    Real Fh[DIMU] = {0.0}; /* flux at neighbour */
    Real R[DIMU][DIMU] = {{0.0}}; /* vector space {Rn} */
    Real Phi[DIMU] = {0.0}; /* flux projection or decomposition coefficients on vector space {Rn} */
    ComputeNonViscousFluxZ(F, k, j, i, U, space, model);
    ComputeNonViscousFluxZ(Fh, k + 1, j, i, U, space, model);
    ComputeEigenvectorSpaceRZ(R, k, j, i, U, space, model);
    ComputeFluxDecompositionCoefficientPhiZ(Phi, k, j, i, U, space, model, dt);
    CalculateReconstructedFlux(Fhat, F, Fh, R, Phi);
    return 0;
}
int TVDFluxY(Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    Real F[DIMU] = {0.0}; /* flux at current node */
    Real Fh[DIMU] = {0.0}; /* flux at neighbour */
    Real R[DIMU][DIMU] = {{0.0}}; /* vector space {Rn} */
    Real Phi[DIMU] = {0.0}; /* flux projection or decomposition coefficients on vector space {Rn} */
    ComputeNonViscousFluxY(F, k, j, i, U, space, model);
    ComputeNonViscousFluxY(Fh, k, j + 1, i, U, space, model);
    ComputeEigenvectorSpaceRY(R, k, j, i, U, space, model);
    ComputeFluxDecompositionCoefficientPhiY(Phi, k, j, i, U, space, model, dt);
    CalculateReconstructedFlux(Fhat, F, Fh, R, Phi);
    return 0;
}
int TVDFluxX(Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    Real F[DIMU] = {0.0}; /* flux at current node */
    Real Fh[DIMU] = {0.0}; /* flux at neighbour */
    Real R[DIMU][DIMU] = {{0.0}}; /* vector space {Rn} */
    Real Phi[DIMU] = {0.0}; /* flux projection or decomposition coefficients on vector space {Rn} */
    ComputeNonViscousFluxX(F, k, j, i, U, space, model);
    ComputeNonViscousFluxX(Fh, k, j, i + 1, U, space, model);
    ComputeEigenvectorSpaceRX(R, k, j, i, U, space, model);
    ComputeFluxDecompositionCoefficientPhiX(Phi, k, j, i, U, space, model, dt);
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
static int ComputeFluxDecompositionCoefficientPhiZ(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    Real g[DIMU] = {0.0}; /* TVD function g at current node */
    Real gh[DIMU] = {0.0}; /* TVD function g at neighbour */
    Real gamma[DIMU] = {0.0}; /* TVD function gamma */
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real alpha[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[DIMU] = {0.0}; /* numerical dissipation */
    ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(lambda, alpha, k, j, i, U, space, model);
    ComputeFunctionGZ(g, k, j, i, U, space, model, dt);
    ComputeFunctionGZ(gh, k + 1, j, i, U, space, model, dt);
    ComputeNumericalDissipationDeltaZ(delta, k, j, i, U, space, model);
    CalculateGamma(gamma, g, gh, alpha);
    for (int row = 0; row < DIMU; ++row) {
        Phi[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta[row]) * alpha[row];
    }
    return 0;
}
static int ComputeFluxDecompositionCoefficientPhiY(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    Real g[DIMU] = {0.0}; /* TVD function g at current node */
    Real gh[DIMU] = {0.0}; /* TVD function g at neighbour */
    Real gamma[DIMU] = {0.0}; /* TVD function gamma */
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real alpha[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[DIMU] = {0.0}; /* numerical dissipation */
    ComputeEigenvaluesAndDecompositionCoefficientAlphaY(lambda, alpha, k, j, i, U, space, model);
    ComputeFunctionGY(g, k, j, i, U, space, model, dt);
    ComputeFunctionGY(gh, k, j + 1, i, U, space, model, dt);
    ComputeNumericalDissipationDeltaY(delta, k, j, i, U, space, model);
    CalculateGamma(gamma, g, gh, alpha);
    for (int row = 0; row < DIMU; ++row) {
        Phi[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta[row]) * alpha[row];
    }
    return 0;
}
static int ComputeFluxDecompositionCoefficientPhiX(
        Real Phi[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    Real g[DIMU] = {0.0}; /* TVD function g at current node */
    Real gh[DIMU] = {0.0}; /* TVD function g at neighbour */
    Real gamma[DIMU] = {0.0}; /* TVD function gamma */
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real alpha[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[DIMU] = {0.0}; /* numerical dissipation */
    ComputeEigenvaluesAndDecompositionCoefficientAlphaX(lambda, alpha, k, j, i, U, space, model);
    ComputeFunctionGX(g, k, j, i, U, space, model, dt);
    ComputeFunctionGX(gh, k, j, i + 1, U, space, model, dt);
    ComputeNumericalDissipationDeltaX(delta, k, j, i, U, space, model);
    CalculateGamma(gamma, g, gh, alpha);
    for (int row = 0; row < DIMU; ++row) {
        Phi[row] = g[row] + gh[row] - Q(lambda[row] + gamma[row], delta[row]) * alpha[row];
    }
    return 0;
}
static int ComputeFunctionGZ(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real lambdah[DIMU] = {0.0}; /* eigenvalues at neighbour */
    Real alpha[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real alphah[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[DIMU] = {0.0}; /* numerical dissipation */
    Real deltah[DIMU] = {0.0}; /* numerical dissipation */
    Real sigma[DIMU] = {0.0}; /* TVD function sigma */
    Real sigmah[DIMU] = {0.0}; /* TVD function sigma at neighbour */
    const Real r = dt * space->ddz;
    ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(lambda, alpha, k, j, i, U, space, model);
    ComputeNumericalDissipationDeltaZ(delta, k, j, i, U, space, model);
    CalculateSigma(sigma, lambda, delta, r);
    ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(lambdah, alphah, k - 1, j, i, U, space, model);
    ComputeNumericalDissipationDeltaZ(deltah, k - 1, j, i, U, space, model);
    CalculateSigma(sigmah, lambdah, deltah, r);
    for (int row = 0; row < DIMU; ++row) {
        g[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
    }
    return 0;
}
static int ComputeFunctionGY(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real lambdah[DIMU] = {0.0}; /* eigenvalues at neighbour */
    Real alpha[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real alphah[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[DIMU] = {0.0}; /* numerical dissipation */
    Real deltah[DIMU] = {0.0}; /* numerical dissipation */
    Real sigma[DIMU] = {0.0}; /* TVD function sigma */
    Real sigmah[DIMU] = {0.0}; /* TVD function sigma at neighbour */
    const Real r = dt * space->ddy;
    ComputeEigenvaluesAndDecompositionCoefficientAlphaY(lambda, alpha, k, j, i, U, space, model);
    ComputeNumericalDissipationDeltaY(delta, k, j, i, U, space, model);
    CalculateSigma(sigma, lambda, delta, r);
    ComputeEigenvaluesAndDecompositionCoefficientAlphaY(lambdah, alphah, k, j - 1, i, U, space, model);
    ComputeNumericalDissipationDeltaY(deltah, k, j - 1, i, U, space, model);
    CalculateSigma(sigmah, lambdah, deltah, r);
    for (int row = 0; row < DIMU; ++row) {
        g[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
    }
    return 0;
}
static int ComputeFunctionGX(
        Real g[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    Real lambda[DIMU] = {0.0}; /* eigenvalues */
    Real lambdah[DIMU] = {0.0}; /* eigenvalues at neighbour */
    Real alpha[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real alphah[DIMU] = {0.0}; /* vector deltaU decomposition coefficients on vector space {Rn} */
    Real delta[DIMU] = {0.0}; /* numerical dissipation */
    Real deltah[DIMU] = {0.0}; /* numerical dissipation */
    Real sigma[DIMU] = {0.0}; /* TVD function sigma */
    Real sigmah[DIMU] = {0.0}; /* TVD function sigma at neighbour */
    const Real r = dt * space->ddx;
    ComputeEigenvaluesAndDecompositionCoefficientAlphaX(lambda, alpha, k, j, i, U, space, model);
    ComputeNumericalDissipationDeltaX(delta, k, j, i, U, space, model);
    CalculateSigma(sigma, lambda, delta, r);
    ComputeEigenvaluesAndDecompositionCoefficientAlphaX(lambdah, alphah, k, j, i - 1, U, space, model);
    ComputeNumericalDissipationDeltaX(deltah, k, j, i - 1, U, space, model);
    CalculateSigma(sigmah, lambdah, deltah, r);
    for (int row = 0; row < DIMU; ++row) {
        g[row] = minmod(sigma[row] * alpha[row], sigmah[row] * alphah[row]);
    }
    return 0;
}
static int CalculateGamma(
        Real gamma[], const Real g[], const Real gh[], const Real alpha[])
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
static int CalculateSigma(
        Real sigma[], const Real lambda[], const Real delta[], const Real r)
{
    for (int row = 0; row < DIMU; ++row) {
        sigma[row] = 0.5 * (Q(lambda[row], delta[row]) - r * lambda[row] * lambda[row]);
    }
    return 0;
}
static int ComputeNumericalDissipationDeltaZ(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    /* numerical dissipation in [0.05, 0.25], 0.125 is recommended */
    ComputeRoeAverageZ(Uo, k, j, i, U, space, model);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    for (int row = 0; row < DIMU; ++row) {
        delta[row] = (fabs(u) + fabs(v) + fabs(w) + c) * model->delta; 
    }
    return 0;
}
static int ComputeNumericalDissipationDeltaY(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    /* numerical dissipation in [0.05, 0.25], 0.125 is recommended */
    ComputeRoeAverageY(Uo, k, j, i, U, space, model);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    for (int row = 0; row < DIMU; ++row) {
        delta[row] = (fabs(u) + fabs(v) + fabs(w) + c) * model->delta; 
    }
    return 0;
}
static int ComputeNumericalDissipationDeltaX(
        Real delta[], const int k, const int j, const int i,
        const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    /* numerical dissipation in [0.05, 0.25], 0.125 is recommended */
    ComputeRoeAverageX(Uo, k, j, i, U, space, model);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    for (int row = 0; row < DIMU; ++row) {
        delta[row] = (fabs(u) + fabs(v) + fabs(w) + c) * model->delta; 
    }
    return 0;
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
    return (sign(x) * max(0, min(fabs(x), y * sign(x))));
}
static int sign(const Real x)
{
    if (0 < x) {
        return 1;
    }
    if (0 > x) {
        return -1;
    }
    return 0;
}
static Real min(const Real x, const Real y)
{
    if (x < y) {
        return x;
    }
    return y;
}
static Real max(const Real x, const Real y)
{
    if (x > y) {
        return x;
    }
    return y;
}
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaZ(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxh = IndexMath(k + 1, j, i, space) * DIMU;;
    Real L[DIMU][DIMU] = {{0.0}}; /* store left eigenvectors */
    const Real deltaU[DIMU] = {
        U[idxh] - U[idx],
        U[idxh+1] - U[idx+1],
        U[idxh+2] - U[idx+2],
        U[idxh+3] - U[idx+3],
        U[idxh+4] - U[idx+4]};
    ComputeEigenvaluesAndEigenvectorSpaceLZ(lambda, L, k, j, i, U, space, model);
    CalculateAlpha(alpha, L, deltaU);
    return 0;
}
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaY(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxh = IndexMath(k, j + 1, i, space) * DIMU;;
    Real L[DIMU][DIMU] = {{0.0}}; /* store left eigenvectors */
    const Real deltaU[DIMU] = {
        U[idxh] - U[idx],
        U[idxh+1] - U[idx+1],
        U[idxh+2] - U[idx+2],
        U[idxh+3] - U[idx+3],
        U[idxh+4] - U[idx+4]};
    ComputeEigenvaluesAndEigenvectorSpaceLY(lambda, L, k, j, i, U, space, model);
    CalculateAlpha(alpha, L, deltaU);
    return 0;
}
static int ComputeEigenvaluesAndDecompositionCoefficientAlphaX(
        Real lambda[], Real alpha[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxh = IndexMath(k, j, i + 1, space) * DIMU;;
    Real L[DIMU][DIMU] = {{0.0}}; /* store left eigenvectors */
    const Real deltaU[DIMU] = {
        U[idxh] - U[idx],
        U[idxh+1] - U[idx+1],
        U[idxh+2] - U[idx+2],
        U[idxh+3] - U[idx+3],
        U[idxh+4] - U[idx+4]};
    ComputeEigenvaluesAndEigenvectorSpaceLX(lambda, L, k, j, i, U, space, model);
    CalculateAlpha(alpha, L, deltaU);
    return 0;
}
static int CalculateAlpha(
        Real alpha[], Real L[][DIMU], const Real deltaU[])
{
    for (int row = 0; row < DIMU; ++row) {
        alpha[row] = 0;
        for (int dummy = 0; dummy < DIMU; ++dummy) {
            alpha[row] = alpha[row] + L[row][dummy] * deltaU[dummy];
        }
    }
    return 0;
}
static int ComputeEigenvaluesAndEigenvectorSpaceLZ(
        Real lambda[], Real L[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageZ(Uo, k, j, i, U, space, model);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (model->gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    lambda[0] = w - c; lambda[1] = w; lambda[2] = w; lambda[3] = w; lambda[4] = w + c;
    L[0][0] = b * q + d * w;   L[0][1] = -b * u;             L[0][2] = -b * v;             L[0][3] = -b * w - d;     L[0][4] = b;
    L[1][0] = -2 * b * q * u;  L[1][1] = 2 * b * u * u + 1;  L[1][2] = 2 * b * v * u;      L[1][3] = 2 * b * w * u;  L[1][4] = -2 * b * u;
    L[2][0] = -2 * b * q * v;  L[2][1] = 2 * b * v * u;      L[2][2] = 2 * b * v * v + 1;  L[2][3] = 2 * b * w * v;  L[2][4] = -2 * b * v;
    L[3][0] = -2 * b * q + 1;  L[3][1] = 2 * b * u;          L[3][2] = 2 * b * v;          L[3][3] = 2 * b * w;      L[3][4] = -2 * b;
    L[4][0] = b * q - d * w;   L[4][1] = -b * u;             L[4][2] = -b * v;             L[4][3] = -b * w + d;     L[4][4] = b;
    return 0;
}
static int ComputeEigenvaluesAndEigenvectorSpaceLY(
        Real lambda[], Real L[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageY(Uo, k, j, i, U, space, model);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (model->gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    lambda[0] = v - c; lambda[1] = v; lambda[2] = v; lambda[3] = v; lambda[4] = v + c;
    L[0][0] = b * q + d * v;    L[0][1] = -b * u;             L[0][2] = -b * v - d;     L[0][3] = -b * w;             L[0][4] = b;
    L[1][0] = -2 * b * q * u;   L[1][1] = 2 * b * u * u + 1;  L[1][2] = 2 * b * v * u;  L[1][3] = 2 * b * w * u;      L[1][4] = -2 * b * u;
    L[2][0] = -2 * b * q + 1;   L[2][1] = 2 * b * u;          L[2][2] = 2 * b * v;      L[2][3] = 2 * b * w;          L[2][4] = -2 * b;
    L[3][0] = -2 * b * q * w;   L[3][1] = 2 * b * w * u;      L[3][2] = 2 * b * w * v;  L[3][3] = 2 * b * w * w + 1;  L[3][4] = -2 * b * w;
    L[4][0] = b * q - d * v;    L[4][1] = -b * u;             L[4][2] = -b * v + d;     L[4][3] = -b * w;             L[4][4] = b;
    return 0;
}
static int ComputeEigenvaluesAndEigenvectorSpaceLX(
        Real lambda[], Real L[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageX(Uo, k, j, i, U, space, model);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (model->gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    lambda[0] = u - c; lambda[1] = u; lambda[2] = u; lambda[3] = u; lambda[4] = u + c;
    L[0][0] = b * q + d * u;   L[0][1] = -b * u - d;     L[0][2] = -b * v;             L[0][3] = -b * w;             L[0][4] = b;
    L[1][0] = -2 * b * q + 1;  L[1][1] = 2 * b * u;      L[1][2] = 2 * b * v;          L[1][3] = 2 * b * w;          L[1][4] = -2 * b;
    L[2][0] = -2 * b * q * v;  L[2][1] = 2 * b * v * u;  L[2][2] = 2 * b * v * v + 1;  L[2][3] = 2 * b * w * v;      L[2][4] = -2 * b * v;
    L[3][0] = -2 * b * q * w;  L[3][1] = 2 * b * w * u;  L[3][2] = 2 * b * w * v;      L[3][3] = 2 * b * w * w + 1;  L[3][4] = -2 * b * w;
    L[4][0] = b * q - d * u;   L[4][1] = -b * u + d;     L[4][2] = -b * v;             L[4][3] = -b * w;             L[4][4] = b;
    return 0;
}
static int ComputeEigenvectorSpaceRZ(
        Real R[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageZ(Uo, k, j, i, U, space, model);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real hT = Uo[4];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    R[0][0] = 1;           R[0][1] = 0;  R[0][2] = 0;  R[0][3] = 1;          R[0][4] = 1;
    R[1][0] = u;           R[1][1] = 1;  R[1][2] = 0;  R[1][3] = 0;          R[1][4] = u;
    R[2][0] = v;           R[2][1] = 0;  R[2][2] = 1;  R[2][3] = 0;          R[2][4] = v;
    R[3][0] = w - c;       R[3][1] = 0;  R[3][2] = 0;  R[3][3] = w;          R[3][4] = w + c;
    R[4][0] = hT - w * c;  R[4][1] = u;  R[4][2] = v;  R[4][3] = w * w - q;  R[4][4] = hT + w * c;
    return 0;
}
static int ComputeEigenvectorSpaceRY(
        Real R[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageY(Uo, k, j, i, U, space, model);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real hT = Uo[4];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    R[0][0] = 1;           R[0][1] = 0;  R[0][2] = 1;          R[0][3] = 0;  R[0][4] = 1;
    R[1][0] = u;           R[1][1] = 1;  R[1][2] = 0;          R[1][3] = 0;  R[1][4] = u;
    R[2][0] = v - c;       R[2][1] = 0;  R[2][2] = v;          R[2][3] = 0;  R[2][4] = v + c;
    R[3][0] = w;           R[3][1] = 0;  R[3][2] = 0;          R[3][3] = 1;  R[3][4] = w;
    R[4][0] = hT - v * c;  R[4][1] = u;  R[4][2] = v * v - q;  R[4][3] = w;  R[4][4] = hT + v * c;
    return 0;
}
static int ComputeEigenvectorSpaceRX(
        Real R[][DIMU], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real Uo[DIMUo] = {0.0}; /* store averaged primitive variables rho, u, v, w, hT, c */
    ComputeRoeAverageX(Uo, k, j, i, U, space, model);
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real hT = Uo[4];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    R[0][0] = 1;           R[0][1] = 1;          R[0][2] = 0;  R[0][3] = 0;  R[0][4] = 1;
    R[1][0] = u - c;       R[1][1] = u;          R[1][2] = 0;  R[1][3] = 0;  R[1][4] = u + c;
    R[2][0] = v;           R[2][1] = 0;          R[2][2] = 1;  R[2][3] = 0;  R[2][4] = v;
    R[3][0] = w;           R[3][1] = 0;          R[3][2] = 0;  R[3][3] = 1;  R[3][4] = w;
    R[4][0] = hT - u * c;  R[4][1] = u * u - q;  R[4][2] = v;  R[4][3] = w;  R[4][4] = hT + u * c;
    return 0;
}
static int ComputeRoeAverageZ(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxh = IndexMath(k + 1, j, i, space) * DIMU;;
    CalculateRoeAverageUo(Uo, idx, idxh, U, model);
    return 0;
}
static int ComputeRoeAverageY(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxh = IndexMath(k, j + 1, i, space) * DIMU;;
    CalculateRoeAverageUo(Uo, idx, idxh, U, model);
    return 0;
}
static int ComputeRoeAverageX(
        Real Uo[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxh = IndexMath(k, j, i + 1, space) * DIMU;;
    CalculateRoeAverageUo(Uo, idx, idxh, U, model);
    return 0;
}
static int CalculateRoeAverageUo(
        Real Uo[], const int idx, const int idxh, 
        const Real *U, const Model *model)
{
    const Real gamma = model->gamma;
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real hT = (U[idx+4] / rho) * gamma - 0.5 * (u * u + v * v + w * w) * (gamma - 1.0);
    const Real rho_h = U[idxh];
    const Real u_h = U[idxh+1] / rho_h;
    const Real v_h = U[idxh+2] / rho_h;
    const Real w_h = U[idxh+3] / rho_h;
    const Real hT_h = (U[idxh+4] / rho_h) * gamma - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h) * (gamma - 1.0);
    const Real D = sqrt(rho_h / rho);
    //Uo[0] = rho * (0.5 * (1.0 + D)) * (0.5 * (1.0 + D)); /* rho average, not required */
    Uo[1] = (u + D * u_h) / (1.0 + D); /* u average */
    Uo[2] = (v + D * v_h) / (1.0 + D); /* v average */
    Uo[3] = (w + D * w_h) / (1.0 + D); /* w average */
    Uo[4] = (hT + D * hT_h) / (1.0 + D); /* hT average */
    Uo[5] = sqrt((gamma - 1.0) * (Uo[4] - 0.5 * (Uo[1] * Uo[1] + Uo[2] * Uo[2] + Uo[3] * Uo[3]))); /* the speed of sound */
    return 0;
}
static int ComputeNonViscousFluxZ(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = rho * (eT - 0.5 * (u * u + v * v + w * w)) * (model->gamma - 1.0);
    F[0] = rho * w;
    F[1] = rho * w * u;
    F[2] = rho * w * v;
    F[3] = rho * w * w + p;
    F[4] = (rho * eT + p) * w;
    return 0;
}
static int ComputeNonViscousFluxY(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = rho * (eT - 0.5 * (u * u + v * v + w * w)) * (model->gamma - 1.0);
    F[0] = rho * v;
    F[1] = rho * v * u;
    F[2] = rho * v * v + p;
    F[3] = rho * v * w;
    F[4] = (rho * eT + p) * v;
    return 0;
}
static int ComputeNonViscousFluxX(
        Real F[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = rho * (eT - 0.5 * (u * u + v * v + w * w)) * (model->gamma - 1.0);
    F[0] = rho * u;
    F[1] = rho * u * u + p;
    F[2] = rho * u * v;
    F[3] = rho * u * w;
    F[4] = (rho * eT + p) * u;
    return 0;
}
/* a good practice: end file with a newline */

