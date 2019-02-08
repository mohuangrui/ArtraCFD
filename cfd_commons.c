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
#include <float.h> /* size of floating point values */
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
typedef void (*EigenvalueSplitter)(const Real [restrict], Real [restrict], Real [restrict]);
typedef void (*EigenvectorLComputer)(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
typedef void (*EigenvectorRComputer)(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
typedef void (*ConvectiveFluxComputer)(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict]);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void LocalLaxFriedrichs(const Real [restrict], Real [restrict], Real [restrict]);
static void StegerWarming(const Real [restrict], Real [restrict], Real [restrict]);
static void EigenvectorLX(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorLY(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorLZ(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorRX(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorRY(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
static void EigenvectorRZ(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict][DIMU]);
static void ConvectiveFluxX(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict]);
static void ConvectiveFluxY(const Real, const Real, const Real, const Real,
        const Real, const Real, Real [restrict]);
static void ConvectiveFluxZ(const Real, const Real, const Real, const Real,
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
    const Real uL = UL[1] / UL[0];
    const Real vL = UL[2] / UL[0];
    const Real wL = UL[3] / UL[0];
    const Real hTL = (UL[4] / UL[0]) * gamma - 0.5 * (uL * uL + vL * vL + wL * wL) * (gamma - 1.0);
    const Real rhoR = UR[0];
    const Real uR = UR[1] / UR[0];
    const Real vR = UR[2] / UR[0];
    const Real wR = UR[3] / UR[0];
    const Real hTR = (UR[4] / UR[0]) * gamma - 0.5 * (uR * uR + vR * vR + wR * wR) * (gamma - 1.0);
    Real D = 0.0;
    switch (averager) {
        case 0: /* arithmetic mean */
            Uo[1] = 0.5 * (uL + uR); /* u average */
            Uo[2] = 0.5 * (vL + vR); /* v average */
            Uo[3] = 0.5 * (wL + wR); /* w average */
            Uo[4] = 0.5 * (hTL + hTR); /* hT average */
            break;
        case 1: /* Roe average */
            D = sqrt(rhoR / rhoL);
            Uo[1] = (uL + D * uR) / (1.0 + D); /* u average */
            Uo[2] = (vL + D * vR) / (1.0 + D); /* v average */
            Uo[3] = (wL + D * wR) / (1.0 + D); /* w average */
            Uo[4] = (hTL + D * hTR) / (1.0 + D); /* hT average */
            break;
        default:
            break;
    }
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
    for (int r = 0; r < DIMU; ++r) {
        LambdaP[r] = 0.5 * (Lambda[r] + lambdaStar);
        LambdaN[r] = 0.5 * (Lambda[r] - lambdaStar);
    }
    return;
}
static void StegerWarming(const Real Lambda[restrict],
        Real LambdaP[restrict], Real LambdaN[restrict])
{
    const Real epsilon = 1.0e-3;
    for (int r = 0; r < DIMU; ++r) {
        LambdaP[r] = 0.5 * (Lambda[r] + sqrt(Lambda[r] * Lambda[r] + epsilon * epsilon));
        LambdaN[r] = 0.5 * (Lambda[r] - sqrt(Lambda[r] * Lambda[r] + epsilon * epsilon));
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
    const Real d = 1.0 / (2.0 * c);
    ComputeEigenvectorL[s](u, v, w, q, b, d, L);
    return;
}
static void EigenvectorLX(const Real u, const Real v, const Real w,
        const Real q, const Real b, const Real d, Real L[restrict][DIMU])
{
    L[0][0] = b * q + d * u;          L[0][1] = -b * u - d;             L[0][2] = -b * v;
    L[0][3] = -b * w;                 L[0][4] = b;
    L[1][0] = -2.0 * b * q + 1.0;     L[1][1] = 2.0 * b * u;            L[1][2] = 2.0 * b * v;
    L[1][3] = 2.0 * b * w;            L[1][4] = -2.0 * b;
    L[2][0] = -2.0 * b * q * v;       L[2][1] = 2.0 * b * v * u;        L[2][2] = 2.0 * b * v * v + 1.0;
    L[2][3] = 2.0 * b * w * v;        L[2][4] = -2.0 * b * v;
    L[3][0] = -2.0 * b * q * w;       L[3][1] = 2.0 * b * w * u;        L[3][2] = 2.0 * b * w * v;
    L[3][3] = 2.0 * b * w * w + 1.0;  L[3][4] = -2.0 * b * w;
    L[4][0] = b * q - d * u;          L[4][1] = -b * u + d;             L[4][2] = -b * v;
    L[4][3] = -b * w;                 L[4][4] = b;
    return;
}
static void EigenvectorLY(const Real u, const Real v, const Real w,
        const Real q, const Real b, const Real d, Real L[restrict][DIMU])
{
    L[0][0] = b * q + d * v;          L[0][1] = -b * u;                 L[0][2] = -b * v - d;
    L[0][3] = -b * w;                 L[0][4] = b;
    L[1][0] = -2.0 * b * q * u;       L[1][1] = 2.0 * b * u * u + 1.0;  L[1][2] = 2.0 * b * v * u;
    L[1][3] = 2.0 * b * w * u;        L[1][4] = -2.0 * b * u;
    L[2][0] = -2.0 * b * q + 1.0;     L[2][1] = 2.0 * b * u;            L[2][2] = 2.0 * b * v;
    L[2][3] = 2.0 * b * w;            L[2][4] = -2.0 * b;
    L[3][0] = -2.0 * b * q * w;       L[3][1] = 2.0 * b * w * u;        L[3][2] = 2.0 * b * w * v;
    L[3][3] = 2.0 * b * w * w + 1.0;  L[3][4] = -2.0 * b * w;
    L[4][0] = b * q - d * v;          L[4][1] = -b * u;                 L[4][2] = -b * v + d;
    L[4][3] = -b * w;                 L[4][4] = b;
    return;
}
static void EigenvectorLZ(const Real u, const Real v, const Real w,
        const Real q, const Real b, const Real d, Real L[restrict][DIMU])
{
    L[0][0] = b * q + d * w;          L[0][1] = -b * u;                 L[0][2] = -b * v;
    L[0][3] = -b * w - d;             L[0][4] = b;
    L[1][0] = -2.0 * b * q * u;       L[1][1] = 2.0 * b * u * u + 1.0;  L[1][2] = 2.0 * b * v * u;
    L[1][3] = 2.0 * b * w * u;        L[1][4] = -2.0 * b * u;
    L[2][0] = -2.0 * b * q * v;       L[2][1] = 2.0 * b * v * u;        L[2][2] = 2.0 * b * v * v + 1.0;
    L[2][3] = 2.0 * b * w * v;        L[2][4] = -2.0 * b * v;
    L[3][0] = -2.0 * b * q + 1.0;     L[3][1] = 2.0 * b * u;            L[3][2] = 2.0 * b * v;
    L[3][3] = 2.0 * b * w;            L[3][4] = -2.0 * b;
    L[4][0] = b * q - d * w;          L[4][1] = -b * u;                 L[4][2] = -b * v;
    L[4][3] = -b * w + d;             L[4][4] = b;
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
static void EigenvectorRX(const Real u, const Real v, const Real w, const Real hT,
        const Real c, const Real q, Real R[restrict][DIMU])
{
    R[0][0] = 1.0;         R[0][1] = 1.0;        R[0][2] = 0.0;  R[0][3] = 0.0;  R[0][4] = 1.0;
    R[1][0] = u - c;       R[1][1] = u;          R[1][2] = 0.0;  R[1][3] = 0.0;  R[1][4] = u + c;
    R[2][0] = v;           R[2][1] = 0.0;        R[2][2] = 1.0;  R[2][3] = 0.0;  R[2][4] = v;
    R[3][0] = w;           R[3][1] = 0.0;        R[3][2] = 0.0;  R[3][3] = 1.0;  R[3][4] = w;
    R[4][0] = hT - u * c;  R[4][1] = u * u - q;  R[4][2] = v;    R[4][3] = w;    R[4][4] = hT + u * c;
    return;
}
static void EigenvectorRY(const Real u, const Real v, const Real w, const Real hT,
        const Real c, const Real q, Real R[restrict][DIMU])
{
    R[0][0] = 1.0;         R[0][1] = 0.0;  R[0][2] = 1.0;        R[0][3] = 0.0;  R[0][4] = 1.0;
    R[1][0] = u;           R[1][1] = 1.0;  R[1][2] = 0.0;        R[1][3] = 0.0;  R[1][4] = u;
    R[2][0] = v - c;       R[2][1] = 0.0;  R[2][2] = v;          R[2][3] = 0.0;  R[2][4] = v + c;
    R[3][0] = w;           R[3][1] = 0.0;  R[3][2] = 0.0;        R[3][3] = 1.0;  R[3][4] = w;
    R[4][0] = hT - v * c;  R[4][1] = u;    R[4][2] = v * v - q;  R[4][3] = w;    R[4][4] = hT + v * c;
    return;
}
static void EigenvectorRZ(const Real u, const Real v, const Real w, const Real hT,
        const Real c, const Real q, Real R[restrict][DIMU])
{
    R[0][0] = 1.0;         R[0][1] = 0.0;  R[0][2] = 0.0;  R[0][3] = 1.0;        R[0][4] = 1.0;
    R[1][0] = u;           R[1][1] = 1.0;  R[1][2] = 0.0;  R[1][3] = 0.0;        R[1][4] = u;
    R[2][0] = v;           R[2][1] = 0.0;  R[2][2] = 1.0;  R[2][3] = 0.0;        R[2][4] = v;
    R[3][0] = w - c;       R[3][1] = 0.0;  R[3][2] = 0.0;  R[3][3] = w;          R[3][4] = w + c;
    R[4][0] = hT - w * c;  R[4][1] = u;    R[4][2] = v;    R[4][3] = w * w - q;  R[4][4] = hT + w * c;
    return;
}
void ConvectiveFlux(const int s, const Real gamma, const Real U[restrict], Real F[restrict])
{
    const Real rho = U[0];
    const Real u = U[1] / U[0];
    const Real v = U[2] / U[0];
    const Real w = U[3] / U[0];
    const Real eT = U[4] / U[0];
    const Real p = ComputePressure(gamma, U);
    ComputeConvectiveFlux[s](rho, u, v, w, eT, p, F);
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
Real Viscosity(const Real T)
{
    return 1.458e-6 * pow(T, 1.5) / (T + 110.4); /* Sutherland's law */
}
Real PrandtlNumber(void)
{
    return 0.71; /* air */
}
void MapPrimitive(const Real gamma, const Real gasR, const Real U[restrict], Real Uo[restrict])
{
    Uo[0] = U[0];
    Uo[1] = U[1] / U[0];
    Uo[2] = U[2] / U[0];
    Uo[3] = U[3] / U[0];
    Uo[4] = (U[4] - 0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / U[0]) * (gamma - 1.0);
    Uo[5] = Uo[4] / (Uo[0] * gasR);
    return;
}
Real ComputePressure(const Real gamma, const Real U[restrict])
{
    return (U[4] - 0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / U[0]) * (gamma - 1.0);
}
Real ComputeTemperature(const Real cv, const Real U[restrict])
{
    return (U[4] - 0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / U[0]) / (U[0] * cv);
}
void MapConservative(const Real gamma, const Real Uo[restrict], Real U[restrict])
{
    U[0] = Uo[0];
    U[1] = Uo[0] * Uo[1];
    U[2] = Uo[0] * Uo[2];
    U[3] = Uo[0] * Uo[3];
    U[4] = 0.5 * Uo[0] * (Uo[1] * Uo[1] + Uo[2] * Uo[2] + Uo[3] * Uo[3]) + Uo[4] / (gamma - 1.0);
    return;
}
int IndexNode(const int k, const int j, const int i, const int jMax, const int iMax)
{
    return (k * jMax + j) * iMax + i;
}
int InPartBox(const int k, const int j, const int i, const int pbox[restrict][LIMIT])
{
    return
        (pbox[Z][MIN] <= k) && (pbox[Z][MAX] > k) &&
        (pbox[Y][MIN] <= j) && (pbox[Y][MAX] > j) &&
        (pbox[X][MIN] <= i) && (pbox[X][MAX] > i);
}
/*
 * Coordinates transformations
 * When transform from spatial coordinates to node coordinates, a half grid
 * distance shift is used to obtain a closest node coordinates, which
 * considers the downward truncation of (int) and works for positive float.
 */
int MapNode(const Real s, const Real sMin, const Real dds, const int ng)
{
    return (int)((s - sMin) * dds + 0.5) + ng;
}
int ConfineSpace(const int n, const int nMin, const int nMax)
{
    return MinInt(nMax - 1, MaxInt(nMin, n));
}
Real MapPoint(const int n, const Real sMin, const Real ds, const int ng)
{
    return sMin + (n - ng) * ds;
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
int EqualReal(const Real x, const Real y)
{
    const Real epsilon = DBL_EPSILON;
    const Real diffMax = FLT_MIN;
    const Real diff = fabs(x - y);
    if (diff <= diffMax) {
        return 1;
    }
    const Real absx = fabs(x);
    const Real absy = fabs(y);
    const Real absMax = (absx > absy) ? absx : absy;
    return (diff <= epsilon * absMax);
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
    const Real zero = 0.0;
    if (zero < x) {
        return 1;
    }
    if (zero > x) {
        return -1;
    }
    return 0;
}
Real Dot(const Real V1[restrict], const Real V2[restrict])
{
    return V1[X] * V2[X] + V1[Y] * V2[Y] + V1[Z] * V2[Z];
}
Real Norm(const Real V[restrict])
{
    return sqrt(Dot(V, V));
}
Real Dist2(const Real V1[restrict], const Real V2[restrict])
{
    const RealVec V = {V2[X] - V1[X], V2[Y] - V1[Y], V2[Z] - V1[Z]};
    return Dot(V, V);
}
Real Dist(const Real V1[restrict], const Real V2[restrict])
{
    return sqrt(Dist2(V1, V2));
}
void Cross(const Real V1[restrict], const Real V2[restrict], Real V[restrict])
{
    V[X] = V1[Y] * V2[Z] - V1[Z] * V2[Y];
    V[Y] = V1[Z] * V2[X] - V1[X] * V2[Z];
    V[Z] = V1[X] * V2[Y] - V1[Y] * V2[X];
    return;
}
void OrthogonalSpace(const Real N[restrict], Real Ta[restrict], Real Tb[restrict])
{
    int mark = Z; /* default mark for minimum component */
    if (fabs(N[mark]) > fabs(N[Y])) {
        mark = Y;
    }
    if (fabs(N[mark]) > fabs(N[X])) {
        mark = X;
    }
    if (X == mark) {
        Ta[X] = 0.0;
        Ta[Y] = -N[Z];
        Ta[Z] = N[Y];
    } else {
        if (Y == mark) {
            Ta[X] = N[Z];
            Ta[Y] = 0.0;
            Ta[Z] = -N[X];
        } else {
            Ta[X] = -N[Y];
            Ta[Y] = N[X];
            Ta[Z] = 0.0;
        }
    }
    Normalize(DIMS, Norm(Ta), Ta);
    Cross(Ta, N, Tb);
    return;
}
void Normalize(const int dimV, const Real normalizer, Real V[restrict])
{
    for (int n = 0; n < dimV; ++n) {
        V[n] = V[n] / normalizer;
    }
    return;
}
/* a good practice: end file with a newline */

