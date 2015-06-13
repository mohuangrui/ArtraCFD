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
    ComputeConvectiveFluxZ(F, k, j, i, U, space, model);
    ComputeConvectiveFluxZ(Fh, k + 1, j, i, U, space, model);
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
    ComputeConvectiveFluxY(F, k, j, i, U, space, model);
    ComputeConvectiveFluxY(Fh, k, j + 1, i, U, space, model);
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
    ComputeConvectiveFluxX(F, k, j, i, U, space, model);
    ComputeConvectiveFluxX(Fh, k, j, i + 1, U, space, model);
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
    return (Sign(x) * MaxReal(0, MinReal(fabs(x), y * Sign(x))));
}
/* a good practice: end file with a newline */

