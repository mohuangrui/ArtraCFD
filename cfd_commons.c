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
#include "cfd_commons.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "commons.h"
/****************************************************************************
 * Function Pointers
 ****************************************************************************/
/*
 * Function pointers are useful for implementing a form of polymorphism.
 * They are mainly used to reduce or avoid switch statement. Pointers to
 * functions can get rather messy. Declaring a typedef to a function pointer
 * generally clarifies the code.
 */
typedef void (*EigenvalueSplitter)(const Real[restrict], Real [restrict], Real [restrict]);
typedef void (*EigenvectorLComputer)(const Real, const Real, const Real, const Real, 
        const Real, const Real, const Real, Real [restrict][DIMU]);
typedef void (*EigenvectorRComputer)(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
typedef void (*ConvectiveFluxComputer)(const Real, const Real, const Real, const Real, 
        const Real, const Real, Real [restrict]);
typedef void (*DiffusiveFluxComputer)(const int, const int, const int, 
        const Space *, const Model *, Real [restrict]);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void LocalLaxFriedrichs(const Real [restrict], Real [restrict], Real [restrict]);
static void StegerWarming(const Real [restrict], Real [restrict], Real [restrict]);
static void EigenvectorLZ(const Real, const Real, const Real, const Real, 
        const Real, const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorLY(const Real, const Real, const Real, const Real, 
        const Real, const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorLX(const Real, const Real, const Real, const Real, 
        const Real, const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorRZ(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorRY(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorRX(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
static void ConvectiveFluxZ(const Real, const Real, const Real, const Real, 
        const Real, const Real, Real [restrict]);
static void ConvectiveFluxY(const Real, const Real, const Real, const Real, 
        const Real, const Real, Real [restrict]);
static void ConvectiveFluxX(const Real, const Real, const Real, const Real, 
        const Real, const Real, Real [restrict]);

/****************************************************************************
 * Global Variables Definition with Private Scope
 ****************************************************************************/
static EigenvalueSplitter SplitEigenvalue[2] = {
    LocalLaxFriedrichs,
    StegerWarming};
static EigenvectorLComputer ComputeEigenvectorL[DIMS] = {
    EigenvectorLX,
    EigenvectorLY,
    EigenvectorLZ};
static EigenvectorRComputer ComputeEigenvectorR[DIMS] = {
    EigenvectorRX,
    EigenvectorRY,
    EigenvectorRZ};
static ConvectiveFluxComputer ComputeConvectiveFlux[DIMS] = {
    ConvectiveFluxX,
    ConvectiveFluxY,
    ConvectiveFluxZ};
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void SymmetricAverage(const int averager, const Real gamma, 
        const Real UL[restrict], const Real UR[restrict], Real Uo[restrict])
{
    const Real rhoL = UL[0];
    const Real uL = UL[1] / rhoL;
    const Real vL = UL[2] / rhoL;
    const Real wL = UL[3] / rhoL;
    const Real hTL = (UL[4] / rhoL) * gamma - 0.5 * (uL * uL + vL * vL + wL * wL) * (gamma - 1.0);
    const Real rhoR = UR[0];
    const Real uR = UR[1] / rhoR;
    const Real vR = UR[2] / rhoR;
    const Real wR = UR[3] / rhoR;
    const Real hTR = (UR[4] / rhoR) * gamma - 0.5 * (uR * uR + vR * vR + wR * wR) * (gamma - 1.0);
    Real D = 1.0; /* default is arithmetic mean */
    if (1 == averager) { /* Roe average */
        D = sqrt(rhoR / rhoL);
    }
    Uo[1] = (uL + D * uR) / (1.0 + D); /* u average */
    Uo[2] = (vL + D * vR) / (1.0 + D); /* v average */
    Uo[3] = (wL + D * wR) / (1.0 + D); /* w average */
    Uo[4] = (hTL + D * hTR) / (1.0 + D); /* hT average */
    Uo[5] = sqrt((gamma - 1.0) * (Uo[4] - 0.5 * (Uo[1] * Uo[1] + Uo[2] * Uo[2] + Uo[3] * Uo[3]))); /* the speed of sound */
    return;
}
void Eigenvalue(const int s, const Real Uo[restrict], Real Lambda[restrict])
{
    Lambda[0] = Uo[s+1] - Uo[5];
    Lambda[1] = Uo[s+1];
    Lambda[2] = Uo[s+1];
    Lambda[3] = Uo[s+1];
    Lambda[4] = Uo[s+1] + Uo[5];
    return;
}
void EigenvalueSplitting(const int splitter, const Real Lambda[restrict],
        Real LambdaP[restrict], Real LambdaN[restrict])
{
    SplitEigenvalue[splitter](Lambda, LambdaP, LambdaN);
    return;
}
static void LocalLaxFriedrichs(const Real Lambda[restrict],
        Real LambdaP[restrict], Real LambdaN[restrict])
{
    /* set local maximum as (|Vs| + c) */
    const Real lambdaStar = fabs(Lambda[2]) + Lambda[4] - Lambda[2];
    for (int row = 0; row < DIMU; ++row) {
        LambdaP[row] = 0.5 * (Lambda[row] + lambdaStar);
        LambdaN[row] = 0.5 * (Lambda[row] - lambdaStar);
    }
    return;
}
static void StegerWarming(const Real Lambda[restrict],
        Real LambdaP[restrict], Real LambdaN[restrict])
{
    const Real epsilon = 1.0e-3;
    for (int row = 0; row < DIMU; ++row) {
        LambdaP[row] = 0.5 * (Lambda[row] + sqrt(Lambda[row] * Lambda[row] + epsilon * epsilon));
        LambdaN[row] = 0.5 * (Lambda[row] - sqrt(Lambda[row] * Lambda[row] + epsilon * epsilon));
    }
    return;
}
void EigenvectorL(const int s, const Real gamma, const Real Uo[restrict], Real L[restrict][DIMU])
{
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    ComputeEigenvectorL[s](u, v, w, c, q, b, d, L);
    return;
}
static void EigenvectorLZ(const Real u, const Real v, const Real w, const Real c, 
        const Real q, const Real b, const Real d, Real L[restrict][DIMU])
{
    L[0][0] = b * q + d * w;   L[0][1] = -b * u;             L[0][2] = -b * v;             L[0][3] = -b * w - d;     L[0][4] = b;
    L[1][0] = -2 * b * q * u;  L[1][1] = 2 * b * u * u + 1;  L[1][2] = 2 * b * v * u;      L[1][3] = 2 * b * w * u;  L[1][4] = -2 * b * u;
    L[2][0] = -2 * b * q * v;  L[2][1] = 2 * b * v * u;      L[2][2] = 2 * b * v * v + 1;  L[2][3] = 2 * b * w * v;  L[2][4] = -2 * b * v;
    L[3][0] = -2 * b * q + 1;  L[3][1] = 2 * b * u;          L[3][2] = 2 * b * v;          L[3][3] = 2 * b * w;      L[3][4] = -2 * b;
    L[4][0] = b * q - d * w;   L[4][1] = -b * u;             L[4][2] = -b * v;             L[4][3] = -b * w + d;     L[4][4] = b;
    return;
}
static void EigenvectorLY(const Real u, const Real v, const Real w, const Real c, 
        const Real q, const Real b, const Real d, Real L[restrict][DIMU])
{
    L[0][0] = b * q + d * v;    L[0][1] = -b * u;             L[0][2] = -b * v - d;     L[0][3] = -b * w;             L[0][4] = b;
    L[1][0] = -2 * b * q * u;   L[1][1] = 2 * b * u * u + 1;  L[1][2] = 2 * b * v * u;  L[1][3] = 2 * b * w * u;      L[1][4] = -2 * b * u;
    L[2][0] = -2 * b * q + 1;   L[2][1] = 2 * b * u;          L[2][2] = 2 * b * v;      L[2][3] = 2 * b * w;          L[2][4] = -2 * b;
    L[3][0] = -2 * b * q * w;   L[3][1] = 2 * b * w * u;      L[3][2] = 2 * b * w * v;  L[3][3] = 2 * b * w * w + 1;  L[3][4] = -2 * b * w;
    L[4][0] = b * q - d * v;    L[4][1] = -b * u;             L[4][2] = -b * v + d;     L[4][3] = -b * w;             L[4][4] = b;
    return;
}
static void EigenvectorLX(const Real u, const Real v, const Real w, const Real c, 
        const Real q, const Real b, const Real d, Real L[restrict][DIMU])
{
    L[0][0] = b * q + d * u;   L[0][1] = -b * u - d;     L[0][2] = -b * v;             L[0][3] = -b * w;             L[0][4] = b;
    L[1][0] = -2 * b * q + 1;  L[1][1] = 2 * b * u;      L[1][2] = 2 * b * v;          L[1][3] = 2 * b * w;          L[1][4] = -2 * b;
    L[2][0] = -2 * b * q * v;  L[2][1] = 2 * b * v * u;  L[2][2] = 2 * b * v * v + 1;  L[2][3] = 2 * b * w * v;      L[2][4] = -2 * b * v;
    L[3][0] = -2 * b * q * w;  L[3][1] = 2 * b * w * u;  L[3][2] = 2 * b * w * v;      L[3][3] = 2 * b * w * w + 1;  L[3][4] = -2 * b * w;
    L[4][0] = b * q - d * u;   L[4][1] = -b * u + d;     L[4][2] = -b * v;             L[4][3] = -b * w;             L[4][4] = b;
    return;
}
void EigenvectorR(const int s, const Real Uo[restrict], Real R[restrict][DIMU])
{
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real hT = Uo[4];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    ComputeEigenvectorR[s](u, v, w, hT, c, q, R);
    return;
}
static void EigenvectorRZ(const Real u, const Real v, const Real w, const Real hT,
        const Real c, const Real q, Real R[restrict][DIMU])
{
    R[0][0] = 1;           R[0][1] = 0;  R[0][2] = 0;  R[0][3] = 1;          R[0][4] = 1;
    R[1][0] = u;           R[1][1] = 1;  R[1][2] = 0;  R[1][3] = 0;          R[1][4] = u;
    R[2][0] = v;           R[2][1] = 0;  R[2][2] = 1;  R[2][3] = 0;          R[2][4] = v;
    R[3][0] = w - c;       R[3][1] = 0;  R[3][2] = 0;  R[3][3] = w;          R[3][4] = w + c;
    R[4][0] = hT - w * c;  R[4][1] = u;  R[4][2] = v;  R[4][3] = w * w - q;  R[4][4] = hT + w * c;
    return;
}
static void EigenvectorRY(const Real u, const Real v, const Real w, const Real hT,
        const Real c, const Real q, Real R[restrict][DIMU])
{
    R[0][0] = 1;           R[0][1] = 0;  R[0][2] = 1;          R[0][3] = 0;  R[0][4] = 1;
    R[1][0] = u;           R[1][1] = 1;  R[1][2] = 0;          R[1][3] = 0;  R[1][4] = u;
    R[2][0] = v - c;       R[2][1] = 0;  R[2][2] = v;          R[2][3] = 0;  R[2][4] = v + c;
    R[3][0] = w;           R[3][1] = 0;  R[3][2] = 0;          R[3][3] = 1;  R[3][4] = w;
    R[4][0] = hT - v * c;  R[4][1] = u;  R[4][2] = v * v - q;  R[4][3] = w;  R[4][4] = hT + v * c;
    return;
}
static void EigenvectorRX(const Real u, const Real v, const Real w, const Real hT,
        const Real c, const Real q, Real R[restrict][DIMU])
{
    R[0][0] = 1;           R[0][1] = 1;          R[0][2] = 0;  R[0][3] = 0;  R[0][4] = 1;
    R[1][0] = u - c;       R[1][1] = u;          R[1][2] = 0;  R[1][3] = 0;  R[1][4] = u + c;
    R[2][0] = v;           R[2][1] = 0;          R[2][2] = 1;  R[2][3] = 0;  R[2][4] = v;
    R[3][0] = w;           R[3][1] = 0;          R[3][2] = 0;  R[3][3] = 1;  R[3][4] = w;
    R[4][0] = hT - u * c;  R[4][1] = u * u - q;  R[4][2] = v;  R[4][3] = w;  R[4][4] = hT + u * c;
    return;
}
void ConvectiveFlux(const int s, const Real gamma, const Real U[restrict], Real F[restrict])
{
    const Real rho = U[0];
    const Real u = U[1] / rho;
    const Real v = U[2] / rho;
    const Real w = U[3] / rho;
    const Real eT = U[4] / rho;
    const Real p = rho * (eT - 0.5 * (u * u + v * v + w * w)) * (gamma - 1.0);
    ComputeConvectiveFlux[s](rho, u, v, w, eT, p, F);
    return;
}
static void ConvectiveFluxZ(const Real rho, const Real u, const Real v, const Real w, 
        const Real eT, const Real p, Real F[restrict])
{
    F[0] = rho * w;
    F[1] = rho * w * u;
    F[2] = rho * w * v;
    F[3] = rho * w * w + p;
    F[4] = (rho * eT + p) * w;
    return;
}
static void ConvectiveFluxY(const Real rho, const Real u, const Real v, const Real w, 
        const Real eT, const Real p, Real F[restrict])
{
    F[0] = rho * v;
    F[1] = rho * v * u;
    F[2] = rho * v * v + p;
    F[3] = rho * v * w;
    F[4] = (rho * eT + p) * v;
    return;
}
static void ConvectiveFluxX(const Real rho, const Real u, const Real v, const Real w, 
        const Real eT, const Real p, Real F[restrict])
{
    F[0] = rho * u;
    F[1] = rho * u * u + p;
    F[2] = rho * u * v;
    F[3] = rho * u * w;
    F[4] = (rho * eT + p) * u;
    return;
}
int DiffusiveFlux(const int s, Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    DiffusiveFluxComputer ComputeDiffusiveFlux[DIMS] = {
        DiffusiveFluxX,
        DiffusiveFluxY,
        DiffusiveFluxZ};
    ComputeDiffusiveFlux[s](G, k, j, i, U, space, model);
    return 0;
}
static void DiffusiveFluxZ(Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxW = IndexMath(k, j, i - 1, space) * DIMU;
    const int idxE = IndexMath(k, j, i + 1, space) * DIMU;
    const int idxS = IndexMath(k, j - 1, i, space) * DIMU;
    const int idxN = IndexMath(k, j + 1, i, space) * DIMU;
    const int idxF = IndexMath(k - 1, j, i, space) * DIMU;
    const int idxB = IndexMath(k + 1, j, i, space) * DIMU;

    /* calculate derivatives in z direction */
    const Real rhoB = U[idxB];
    const Real uB = U[idxB+1] / rhoB;
    const Real vB = U[idxB+2] / rhoB;
    const Real wB = U[idxB+3] / rhoB;
    const Real eTB = U[idxB+4] / rhoB;
    const Real TB = (eTB - 0.5 * (uB * uB + vB * vB + wB * wB)) / model->cv;

    const Real rhoF = U[idxF];
    const Real uF = U[idxF+1] / rhoF;
    const Real vF = U[idxF+2] / rhoF;
    const Real wF = U[idxF+3] / rhoF;
    const Real eTF = U[idxF+4] / rhoF;
    const Real TF = (eTF - 0.5 * (uF * uF + vF * vF + wF * wF)) / model->cv;

    const Real du_dz = (uB - uF) * (0.5 * space->ddz);
    const Real dv_dz = (vB - vF) * (0.5 * space->ddz);
    const Real dw_dz = (wB - wF) * (0.5 * space->ddz);
    const Real dT_dz = (TB - TF) * (0.5 * space->ddz);

    /* calculate derivatives in y direction */
    const Real vN = U[idxN+2] / U[idxN];
    const Real wN = U[idxN+3] / U[idxN];
    const Real vS = U[idxS+2] / U[idxS];
    const Real wS = U[idxS+3] / U[idxS];
    const Real dv_dy = (vN - vS) * (0.5 * space->ddy);
    const Real dw_dy = (wN - wS) * (0.5 * space->ddy);

    /* calculate derivatives in x direction */
    const Real uE = U[idxE+1] / U[idxE];
    const Real wE = U[idxE+3] / U[idxE];
    const Real uW = U[idxW+1] / U[idxW];
    const Real wW = U[idxW+3] / U[idxW];
    const Real du_dx = (uE - uW) * (0.5 * space->ddx);
    const Real dw_dx = (wE - wW) * (0.5 * space->ddx);

    /* the primitive variables in current point */
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real T = (eT - 0.5 * (u * u + v * v + w * w)) / model->cv;

    /* Calculate dynamic viscosity and heat conductivity */
    const Real mu = model->refMu * 1.45e-6 * (pow(T * model->refTemperature, 1.5) / (T * model->refTemperature + 110));
    const Real heatK = model->gamma * model->cv * mu / model->refPr;
    const Real divV = du_dx + dv_dy + dw_dz;

    G[0] = 0;
    G[1] = mu * (dw_dx + du_dz);
    G[2] = mu * (dw_dy + dv_dz);
    G[3] = mu * (2.0 * dw_dz - (2.0/3.0) * divV);
    G[4] = heatK * dT_dz + u * G[1] + v * G[2] + w * G[3];
    return 0;
}
static void DiffusiveFluxY(Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxW = IndexMath(k, j, i - 1, space) * DIMU;
    const int idxE = IndexMath(k, j, i + 1, space) * DIMU;
    const int idxS = IndexMath(k, j - 1, i, space) * DIMU;
    const int idxN = IndexMath(k, j + 1, i, space) * DIMU;
    const int idxF = IndexMath(k - 1, j, i, space) * DIMU;
    const int idxB = IndexMath(k + 1, j, i, space) * DIMU;

    /* calculate derivatives in z direction */
    const Real vB = U[idxB+2] / U[idxB];
    const Real wB = U[idxB+3] / U[idxB];
    const Real vF = U[idxF+2] / U[idxF];
    const Real wF = U[idxF+3] / U[idxF];
    const Real dv_dz = (vB - vF) * (0.5 * space->ddz);
    const Real dw_dz = (wB - wF) * (0.5 * space->ddz);

    /* calculate derivatives in y direction */
    const Real rhoN = U[idxN];
    const Real uN = U[idxN+1] / rhoN;
    const Real vN = U[idxN+2] / rhoN;
    const Real wN = U[idxN+3] / rhoN;
    const Real eTN = U[idxN+4] / rhoN;
    const Real TN = (eTN - 0.5 * (uN * uN + vN * vN + wN * wN)) / model->cv;

    const Real rhoS = U[idxS];
    const Real uS = U[idxS+1] / rhoS;
    const Real vS = U[idxS+2] / rhoS;
    const Real wS = U[idxS+3] / rhoS;
    const Real eTS = U[idxS+4] / rhoS;
    const Real TS = (eTS - 0.5 * (uS * uS + vS * vS + wS * wS)) / model->cv;

    const Real du_dy = (uN - uS) * (0.5 * space->ddy);
    const Real dv_dy = (vN - vS) * (0.5 * space->ddy);
    const Real dw_dy = (wN - wS) * (0.5 * space->ddy);
    const Real dT_dy = (TN - TS) * (0.5 * space->ddy);

    /* calculate derivatives in x direction */
    const Real uE = U[idxE+1] / U[idxE];
    const Real vE = U[idxE+2] / U[idxE];
    const Real uW = U[idxW+1] / U[idxW];
    const Real vW = U[idxW+2] / U[idxW];
    const Real du_dx = (uE - uW) * (0.5 * space->ddx);
    const Real dv_dx = (vE - vW) * (0.5 * space->ddx);

    /* the primitive variables in current point */
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real T = (eT - 0.5 * (u * u + v * v + w * w)) / model->cv;

    /* Calculate dynamic viscosity and heat conductivity */
    const Real mu = model->refMu * 1.45e-6 * (pow(T * model->refTemperature, 1.5) / (T * model->refTemperature + 110));
    const Real heatK = model->gamma * model->cv * mu / model->refPr;
    const Real divV = du_dx + dv_dy + dw_dz;

    G[0] = 0;
    G[1] = mu * (dv_dx + du_dy);
    G[2] = mu * (2.0 * dv_dy - (2.0/3.0) * divV);
    G[3] = mu * (dv_dz + dw_dy);
    G[4] = heatK * dT_dy + u * G[1] + v * G[2] + w * G[3];
    return 0;
}
static void DiffusiveFluxX(Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    const int idx = IndexMath(k, j, i, space) * DIMU;
    const int idxW = IndexMath(k, j, i - 1, space) * DIMU;
    const int idxE = IndexMath(k, j, i + 1, space) * DIMU;
    const int idxS = IndexMath(k, j - 1, i, space) * DIMU;
    const int idxN = IndexMath(k, j + 1, i, space) * DIMU;
    const int idxF = IndexMath(k - 1, j, i, space) * DIMU;
    const int idxB = IndexMath(k + 1, j, i, space) * DIMU;

    /* calculate derivatives in z direction */
    const Real uB = U[idxB+1] / U[idxB];
    const Real uF = U[idxF+1] / U[idxF];
    const Real wB = U[idxB+3] / U[idxB];
    const Real wF = U[idxF+3] / U[idxF];
    const Real du_dz = (uB - uF) * (0.5 * space->ddz);
    const Real dw_dz = (wB - wF) * (0.5 * space->ddz);

    /* calculate derivatives in y direction */
    const Real uN = U[idxN+1] / U[idxN];
    const Real uS = U[idxS+1] / U[idxS];
    const Real vN = U[idxN+2] / U[idxN];
    const Real vS = U[idxS+2] / U[idxS];
    const Real du_dy = (uN - uS) * (0.5 * space->ddy);
    const Real dv_dy = (vN - vS) * (0.5 * space->ddy);

    /* calculate derivatives in x direction */
    const Real rhoE = U[idxE];
    const Real uE = U[idxE+1] / rhoE;
    const Real vE = U[idxE+2] / rhoE;
    const Real wE = U[idxE+3] / rhoE;
    const Real eTE = U[idxE+4] / rhoE;
    const Real TE = (eTE - 0.5 * (uE * uE + vE * vE + wE * wE)) / model->cv;

    const Real rhoW = U[idxW];
    const Real uW = U[idxW+1] / rhoW;
    const Real vW = U[idxW+2] / rhoW;
    const Real wW = U[idxW+3] / rhoW;
    const Real eTW = U[idxW+4] / rhoW;
    const Real TW = (eTW - 0.5 * (uW * uW + vW * vW + wW * wW)) / model->cv;

    const Real du_dx = (uE - uW) * (0.5 * space->ddx);
    const Real dv_dx = (vE - vW) * (0.5 * space->ddx);
    const Real dw_dx = (wE - wW) * (0.5 * space->ddx);
    const Real dT_dx = (TE - TW) * (0.5 * space->ddx);

    /* the primitive variables in current point */
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real T = (eT - 0.5 * (u * u + v * v + w * w)) / model->cv;

    /* Calculate dynamic viscosity and heat conductivity */
    const Real mu = model->refMu * 1.45e-6 * (pow(T * model->refTemperature, 1.5) / (T * model->refTemperature + 110));
    const Real heatK = model->gamma * model->cv * mu / model->refPr;
    const Real divV = du_dx + dv_dy + dw_dz;

    G[0] = 0;
    G[1] = mu * (2.0 * du_dx - (2.0/3.0) * divV);
    G[2] = mu * (du_dy + dv_dx);
    G[3] = mu * (du_dz + dw_dx);
    G[4] = heatK * dT_dx + u * G[1] + v * G[2] + w * G[3];
    return 0;
}
/*
 * Get value of primitive variable vector.
 */
int PrimitiveByConservative(Real Uo[], const int idx, const Real *U, const Model *model)
{
    Uo[0] = U[idx];
    Uo[1] = U[idx+1] / U[idx];
    Uo[2] = U[idx+2] / U[idx];
    Uo[3] = U[idx+3] / U[idx];
    Uo[4] = (U[idx+4] - 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx]) * (model->gamma - 1.0);
    Uo[5] = Uo[4] / (Uo[0] * model->gasR);
    return 0;
}
/*
 * Compute conservative variable vector according to primitives.
 */
int ConservativeByPrimitive(Real *U, const int idx, const Real Uo[], const Model *model)
{
    U[idx] = Uo[0];
    U[idx+1] = Uo[0] * Uo[1];
    U[idx+2] = Uo[0] * Uo[2];
    U[idx+3] = Uo[0] * Uo[3];
    U[idx+4] = 0.5 * Uo[0] * (Uo[1] * Uo[1] + Uo[2] * Uo[2] + Uo[3] * Uo[3]) + Uo[4] / (model->gamma - 1.0); 
    return 0;
}
Real ComputePressure(const int idx, const Real *U, const Model *model)
{
    return (U[idx+4] - 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx]) * (model->gamma - 1.0);
}
Real ComputeTemperature(const int idx, const Real *U, const Model *model)
{
    return (U[idx+4] - 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx]) / (U[idx] * model->cv);
}
/*
 * Index math.
 */
int IndexNode(const int k, const int j, const int i, const Space *space)
{
    return ((k * space->n[Y] + j) * space->n[X] + i);
}
/*
 * Coordinates transformations.
 * When transform from spatial coordinates to node coordinates, a half
 * grid distance shift is necessary to ensure obtaining a closest node
 * coordinates considering the downward truncation of (int). 
 * Note: current rounding conversion only works for positive float.
 */
int NodeSpace(const Real point, const int s, const Space *space)
{
    return (int)((point - space->domain[s][MIN]) * space->dd[s] + 0.5) + space->ng;
}
int ValidNodeSpace(const int node, const int s, const Partition *part)
{
    return MinInt(part->n[PIN][s][MAX] - 1, MaxInt(part->n[PIN][s][MIN], node));
}
Real PointSpace(const int node, const int s, const Space *space)
{
    return (space->domain[s][MIN] + (node - space->ng) * space->d[s]);
}
/*
 * Math functions
 */
Real MinReal(const Real x, const Real y)
{
    if (x < y) {
        return x;
    }
    return y;
}
Real MaxReal(const Real x, const Real y)
{
    if (x > y) {
        return x;
    }
    return y;
}
int MinInt(const int x, const int y)
{
    if (x < y) {
        return x;
    }
    return y;
}
int MaxInt(const int x, const int y)
{
    if (x > y) {
        return x;
    }
    return y;
}
int Sign(const Real x)
{
    if (0 < x) {
        return 1;
    }
    if (0 > x) {
        return -1;
    }
    return 0;
}
Real Dot(const RealVector V1, const RealVector V2)
{
    return (V1[X] * V2[X] + V1[Y] * V2[Y] + V1[Z] * V2[Z]);
}
Real Norm(const RealVector V)
{
    return sqrt(Dot(V, V));
}
Real Dist2(const RealVector V1, const RealVector V2)
{
    const RealVector V = {V1[X] - V2[X], V1[Y] - V2[Y], V1[Z] - V2[Z]};
    return Dot(V, V);
}
Real Dist(const RealVector V1, const RealVector V2)
{
    return sqrt(Dist2(V1, V2));
}
int Cross(RealVector V, const RealVector V1, const RealVector V2)
{
    V[X] = V1[Y] * V2[Z] - V1[Z] * V2[Y];
    V[Y] = V1[Z] * V2[X] - V1[X] * V2[Z];
    V[Z] = V1[X] * V2[Y] - V1[Y] * V2[X];
    return 0;
}
int OrthogonalSpace(const RealVector N, RealVector Ta, RealVector Tb)
{
    int mark = Z; /* default mark for minimum component */
    if (fabs(N[mark]) > fabs(N[Y])) {
        mark = Y;
    }
    if (fabs(N[mark]) > fabs(N[X])) {
        mark = X;
    }
    if (X == mark) {
        Ta[X] = 0;
        Ta[Y] = -N[Z];
        Ta[Z] = N[Y];
    } else {
        if (Y == mark) {
            Ta[X] = -N[Z];
            Ta[Y] = 0;
            Ta[Z] = N[X];
        } else {
            Ta[X] = -N[Y];
            Ta[Y] = N[X];
            Ta[Z] = 0;
        }
    }
    Normalize(Ta, DIMS, Norm(Ta));
    Cross(Tb, N, Ta);
    return 0;
}
int Normalize(Real V[], const int dimV, const Real normalizer)
{
    for (int n = 0; n < dimV; ++n) {
        V[n] = V[n] / normalizer;
    }
    return 0;
}
/* a good practice: end file with a newline */

