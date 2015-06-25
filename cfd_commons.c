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
 * Function Pointers
 ****************************************************************************/
/*
 * Function pointers are useful for implementing a form of polymorphism.
 * They are mainly used to reduce or avoid switch statement. Pointers to
 * functions can get rather messy. Declaring a typedef to a function pointer
 * generally clarifies the code.
 */
typedef int (*EigenvectorSpaceLComputer)(Real [][DIMU], const Real [], const Real);
typedef int (*EigenvectorSpaceRComputer)(Real [][DIMU], const Real []);
typedef int (*ConvectiveFluxComputer)(Real [], const int, const Real *, const Real);
typedef int (*DiffusiveFluxComputer)(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int EigenvectorSpaceLZ(Real [][DIMU], const Real [], const Real);
static int EigenvectorSpaceLY(Real [][DIMU], const Real [], const Real);
static int EigenvectorSpaceLX(Real [][DIMU], const Real [], const Real);
static int EigenvectorSpaceRZ(Real [][DIMU], const Real []);
static int EigenvectorSpaceRY(Real [][DIMU], const Real []);
static int EigenvectorSpaceRX(Real [][DIMU], const Real []);
static int ConvectiveFluxZ(Real [], const int, const Real *, const Real);
static int ConvectiveFluxY(Real [], const int, const Real *, const Real);
static int ConvectiveFluxX(Real [], const int, const Real *, const Real);
static int DiffusiveFluxZ(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
static int DiffusiveFluxY(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
static int DiffusiveFluxX(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int DecompositionCoefficientAlpha(Real alpha[], Real L[][DIMU], const int idx, const int idxh, const Real *U)
{
    const Real deltaU[DIMU] = { /* U variation */
        U[idxh] - U[idx],
        U[idxh+1] - U[idx+1],
        U[idxh+2] - U[idx+2],
        U[idxh+3] - U[idx+3],
        U[idxh+4] - U[idx+4]};
    for (int row = 0; row < DIMU; ++row) {
        alpha[row] = 0;
        for (int dummy = 0; dummy < DIMU; ++dummy) {
            alpha[row] = alpha[row] + L[row][dummy] * deltaU[dummy];
        }
    }
    return 0;
}
int EigenvalueLambda(const int s, Real lambda[], const Real Uo[])
{
    lambda[0] = Uo[s+1] - Uo[5];
    lambda[1] = Uo[s+1];
    lambda[2] = Uo[s+1];
    lambda[3] = Uo[s+1];
    lambda[4] = Uo[s+1] + Uo[5];
    return 0;
}
int EigenvectorSpaceL(const int s, Real L[][DIMU], const Real Uo[], const Real gamma)
{
    EigenvectorSpaceLComputer ComputeEigenvectorSpaceL[DIMS] = {
        EigenvectorSpaceLX,
        EigenvectorSpaceLY,
        EigenvectorSpaceLZ};
    return 0;
}
static int EigenvectorSpaceLZ(Real L[][DIMU], const Real Uo[], const Real gamma)
{
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    L[0][0] = b * q + d * w;   L[0][1] = -b * u;             L[0][2] = -b * v;             L[0][3] = -b * w - d;     L[0][4] = b;
    L[1][0] = -2 * b * q * u;  L[1][1] = 2 * b * u * u + 1;  L[1][2] = 2 * b * v * u;      L[1][3] = 2 * b * w * u;  L[1][4] = -2 * b * u;
    L[2][0] = -2 * b * q * v;  L[2][1] = 2 * b * v * u;      L[2][2] = 2 * b * v * v + 1;  L[2][3] = 2 * b * w * v;  L[2][4] = -2 * b * v;
    L[3][0] = -2 * b * q + 1;  L[3][1] = 2 * b * u;          L[3][2] = 2 * b * v;          L[3][3] = 2 * b * w;      L[3][4] = -2 * b;
    L[4][0] = b * q - d * w;   L[4][1] = -b * u;             L[4][2] = -b * v;             L[4][3] = -b * w + d;     L[4][4] = b;
    return 0;
}
static int EigenvectorSpaceLY(Real L[][DIMU], const Real Uo[], const Real gamma)
{
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    L[0][0] = b * q + d * v;    L[0][1] = -b * u;             L[0][2] = -b * v - d;     L[0][3] = -b * w;             L[0][4] = b;
    L[1][0] = -2 * b * q * u;   L[1][1] = 2 * b * u * u + 1;  L[1][2] = 2 * b * v * u;  L[1][3] = 2 * b * w * u;      L[1][4] = -2 * b * u;
    L[2][0] = -2 * b * q + 1;   L[2][1] = 2 * b * u;          L[2][2] = 2 * b * v;      L[2][3] = 2 * b * w;          L[2][4] = -2 * b;
    L[3][0] = -2 * b * q * w;   L[3][1] = 2 * b * w * u;      L[3][2] = 2 * b * w * v;  L[3][3] = 2 * b * w * w + 1;  L[3][4] = -2 * b * w;
    L[4][0] = b * q - d * v;    L[4][1] = -b * u;             L[4][2] = -b * v + d;     L[4][3] = -b * w;             L[4][4] = b;
    return 0;
}
static int EigenvectorSpaceLX(Real L[][DIMU], const Real Uo[], const Real gamma)
{
    const Real u = Uo[1];
    const Real v = Uo[2];
    const Real w = Uo[3];
    const Real c = Uo[5];
    const Real q = 0.5 * (u * u + v * v + w * w);
    const Real b = (gamma - 1.0) / (2.0 * c * c);
    const Real d = (1.0 / (2.0 * c)); 
    L[0][0] = b * q + d * u;   L[0][1] = -b * u - d;     L[0][2] = -b * v;             L[0][3] = -b * w;             L[0][4] = b;
    L[1][0] = -2 * b * q + 1;  L[1][1] = 2 * b * u;      L[1][2] = 2 * b * v;          L[1][3] = 2 * b * w;          L[1][4] = -2 * b;
    L[2][0] = -2 * b * q * v;  L[2][1] = 2 * b * v * u;  L[2][2] = 2 * b * v * v + 1;  L[2][3] = 2 * b * w * v;      L[2][4] = -2 * b * v;
    L[3][0] = -2 * b * q * w;  L[3][1] = 2 * b * w * u;  L[3][2] = 2 * b * w * v;      L[3][3] = 2 * b * w * w + 1;  L[3][4] = -2 * b * w;
    L[4][0] = b * q - d * u;   L[4][1] = -b * u + d;     L[4][2] = -b * v;             L[4][3] = -b * w;             L[4][4] = b;
    return 0;
}
int EigenvectorSpaceR(const int s, Real R[][DIMU], const Real Uo[])
{
    EigenvectorSpaceRComputer ComputeEigenvectorSpaceR[DIMS] = {
        EigenvectorSpaceRX,
        EigenvectorSpaceRY,
        EigenvectorSpaceRZ};
    return 0;
}
static int EigenvectorSpaceRZ(Real R[][DIMU], const Real Uo[])
{
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
static int EigenvectorSpaceRY(Real R[][DIMU], const Real Uo[])
{
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
static int EigenvectorSpaceRX(Real R[][DIMU], const Real Uo[])
{
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
int ComputeRoeAverage(Real Uo[], const int idx, const int idxh, const Real *U, const Real gamma)
{
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
int ConvectiveFlux(const int s, Real F[], const int idx, const Real *U, const Real gamma)
{
    ConvectiveFluxComputer ComputeConvectiveFlux[DIMS] = {
        ConvectiveFluxX,
        ConvectiveFluxY,
        ConvectiveFluxZ};
    return 0;
}
static int ConvectiveFluxZ(Real F[], const int idx, const Real *U, const Real gamma)
{
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = rho * (eT - 0.5 * (u * u + v * v + w * w)) * (gamma - 1.0);
    F[0] = rho * w;
    F[1] = rho * w * u;
    F[2] = rho * w * v;
    F[3] = rho * w * w + p;
    F[4] = (rho * eT + p) * w;
    return 0;
}
static int ConvectiveFluxY(Real F[], const int idx, const Real *U, const Real gamma)
{
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = rho * (eT - 0.5 * (u * u + v * v + w * w)) * (gamma - 1.0);
    F[0] = rho * v;
    F[1] = rho * v * u;
    F[2] = rho * v * v + p;
    F[3] = rho * v * w;
    F[4] = (rho * eT + p) * v;
    return 0;
}
static int ConvectiveFluxX(Real F[], const int idx, const Real *U, const Real gamma)
{
    const Real rho = U[idx];
    const Real u = U[idx+1] / rho;
    const Real v = U[idx+2] / rho;
    const Real w = U[idx+3] / rho;
    const Real eT = U[idx+4] / rho;
    const Real p = rho * (eT - 0.5 * (u * u + v * v + w * w)) * (gamma - 1.0);
    F[0] = rho * u;
    F[1] = rho * u * u + p;
    F[2] = rho * u * v;
    F[3] = rho * u * w;
    F[4] = (rho * eT + p) * u;
    return 0;
}
int DiffusiveFlux(const int s, Real G[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    DiffusiveFluxComputer ComputeDiffusiveFlux[DIMS] = {
        DiffusiveFluxX,
        DiffusiveFluxY,
        DiffusiveFluxZ};
    return 0;
}
static int DiffusiveFluxZ(Real G[], const int k, const int j, const int i, 
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
static int DiffusiveFluxY(Real G[], const int k, const int j, const int i, 
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
static int DiffusiveFluxX(Real G[], const int k, const int j, const int i, 
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
int IndexMath(const int k, const int j, const int i, const Space *space)
{
    return ((k * space->jMax + j) * space->iMax + i);
}
Real *IndexGeometry(const int geoID, const Geometry *geometry)
{
    return geometry->headAddress + geoID * ENTRYGEO;
}
/*
 * Coordinates transformations.
 * When transform from spatial coordinates to node coordinates, a half
 * grid distance shift is necessary to ensure obtaining a closest node
 * coordinates considering the downward truncation of (int). 
 * Note: current rounding conversion only works for positive float.
 */
int ComputeK(const Real z, const Space *space)
{
    return (int)((z - space->zMin) * space->ddz + 0.5) + space->ng;
}
int ComputeJ(const Real y, const Space *space)
{
    return (int)((y - space->yMin) * space->ddy + 0.5) + space->ng;
}
int ComputeI(const Real x, const Space *space)
{
    return (int)((x - space->xMin) * space->ddx + 0.5) + space->ng;
}
int ValidRegionK(const int k, const Partition *part)
{
    return MinInt(part->kSup[0] - 1, MaxInt(part->kSub[0], k));
}
int ValidRegionJ(const int j, const Partition *part)
{
    return MinInt(part->jSup[0] - 1, MaxInt(part->jSub[0], j));
}
int ValidRegionI(const int i, const Partition *part)
{
    return MinInt(part->iSup[0] - 1, MaxInt(part->iSub[0], i));
}
Real ComputeZ(const int k, const Space *space)
{
    return (space->zMin + (k - space->ng) * space->dz);
}
Real ComputeY(const int j, const Space *space)
{
    return (space->yMin + (j - space->ng) * space->dy);
}
Real ComputeX(const int i, const Space *space)
{
    return (space->xMin + (i - space->ng) * space->dx);
}
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
/* a good practice: end file with a newline */

