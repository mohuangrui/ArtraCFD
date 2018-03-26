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
#include "fluid_dynamics.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "weno.h"
#include "boundary_treatment.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Function Pointers
 ****************************************************************************/
typedef void (*TemporalDiscretizer)(const Real, const int, Space *, const Model *);
typedef void (*SolutionOperator)(const int, const Real, const Real, 
        const Real [restrict], const Real [restrict], Real [restrict], const Real [restrict], 
        const Real [restrict], const Real [restrict], const Real [restrict],
        const Real [restrict], const Real, const Real [restrict]);
typedef void (*ConvectiveFluxReconstructor)(const int, const int, const int,
        const int, const int, const int [restrict], const Node *const, 
        const Model *, Real [restrict]);
typedef void (*DiffusiveFluxReconstructor)(const int, const int, const int, 
        const int, const int [restrict], const Real [restrict], const Node *const, 
        const Model *, Real [restrict]);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void TemporalDiscretization(const Real, const int, Space *, const Model *);
static void RungeKutta2(const Real, const int, Space *, const Model *);
static void RungeKutta3(const Real, const int, Space *, const Model *);
static void LLL(const Real, const Real, const Real, const int,
        const int, const int, const int, Space *, const Model *);
static void DimensionalSplitting(const int, const Real, const Real, 
        const Real [restrict], const Real [restrict], Real [restrict], const Real [restrict], 
        const Real [restrict], const Real [restrict], const Real [restrict],
        const Real [restrict], const Real, const Real [restrict]);
static void DimensionByDimension(const int, const Real, const Real, 
        const Real [restrict], const Real [restrict], Real [restrict], const Real [restrict], 
        const Real [restrict], const Real [restrict], const Real [restrict],
        const Real [restrict], const Real, const Real [restrict]);
static void NumericalConvectiveFlux(const int, const int, const int, const int,
        const int, const int [restrict], const Node *const, const Model *, Real [restrict]);
static void NumericalDiffusiveFlux(const int, const int, const int, const int, 
        const int, const int [restrict], const Real [restrict], const Node *const, 
        const Model *, Real [restrict]);
static void NumericalDiffusiveFluxX(const int, const int, const int, const int,
        const int [restrict], const Real [restrict], const Node *const, 
        const Model *, Real [restrict]);
static void NumericalDiffusiveFluxY(const int, const int, const int, const int,
        const int [restrict], const Real [restrict], const Node *const, 
        const Model *, Real [restrict]);
static void NumericalDiffusiveFluxZ(const int, const int, const int, const int,
        const int [restrict], const Real [restrict], const Node *const, 
        const Model *, Real [restrict]);
static void SourceVector(const int, const int, const int, const int,
        const int [restrict], const Node *const, const Model *, Real [restrict]);
/****************************************************************************
 * Global Variables Definition with Private Scope
 ****************************************************************************/
static TemporalDiscretizer DiscretizeTime[2] = {
    RungeKutta2,
    RungeKutta3};
static SolutionOperator SolveOperator[2] = {
    DimensionalSplitting,
    DimensionByDimension};
static ConvectiveFluxReconstructor ReconstructConvectiveFlux[2] = {
    WENO3,
    WENO5};
static DiffusiveFluxReconstructor ReconstructDiffusiveFlux[DIMS] = {
    NumericalDiffusiveFluxX,
    NumericalDiffusiveFluxY,
    NumericalDiffusiveFluxZ};
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Two approaches for addressing multidimensional governing equations:
 *   a) - dimensional splitting
 *   b) - dimension-by-dimension approximation
 * Time and space discretizations are implemented under the method of lines 
 */
void FluidDynamics(const Real dt, Space *space, const Model *model)
{
    switch (model->multidim) {
        case 0: /* dimensional splitting approximation */
            break;
        case 1: /* dimension-by-dimension approximation */
            TemporalDiscretization(0.5 * dt, DIMS, space, model);
            TemporalDiscretization(0.5 * dt, DIMS, space, model);
            return;
        default:
            break;
    }
    switch (space->part.collapse) {
        case COLLAPSEN:
            TemporalDiscretization(0.5 * dt, Z, space, model);
            TemporalDiscretization(0.5 * dt, Y, space, model);
            TemporalDiscretization(0.5 * dt, X, space, model);
            TemporalDiscretization(0.5 * dt, X, space, model);
            TemporalDiscretization(0.5 * dt, Y, space, model);
            TemporalDiscretization(0.5 * dt, Z, space, model);
            break;
        case COLLAPSEX:
            TemporalDiscretization(0.5 * dt, Z, space, model);
            TemporalDiscretization(0.5 * dt, Y, space, model);
            TemporalDiscretization(0.5 * dt, Y, space, model);
            TemporalDiscretization(0.5 * dt, Z, space, model);
            break;
        case COLLAPSEY:
            TemporalDiscretization(0.5 * dt, Z, space, model);
            TemporalDiscretization(0.5 * dt, X, space, model);
            TemporalDiscretization(0.5 * dt, X, space, model);
            TemporalDiscretization(0.5 * dt, Z, space, model);
            break;
        case COLLAPSEZ:
            TemporalDiscretization(0.5 * dt, Y, space, model);
            TemporalDiscretization(0.5 * dt, X, space, model);
            TemporalDiscretization(0.5 * dt, X, space, model);
            TemporalDiscretization(0.5 * dt, Y, space, model);
            break;
        case COLLAPSEXY:
            TemporalDiscretization(0.5 * dt, Z, space, model);
            TemporalDiscretization(0.5 * dt, Z, space, model);
            break;
        case COLLAPSEXZ:
            TemporalDiscretization(0.5 * dt, Y, space, model);
            TemporalDiscretization(0.5 * dt, Y, space, model);
            break;
        case COLLAPSEYZ:
            TemporalDiscretization(0.5 * dt, X, space, model);
            TemporalDiscretization(0.5 * dt, X, space, model);
            break;
        default:
            break;
    }
    return;
}
/*
 * Computation must start from TO data space and end with TO data space.
 */
static void TemporalDiscretization(const Real dt, const int s, Space *space, const Model *model)
{
    DiscretizeTime[model->tScheme](dt, s, space, model);
    return;
}
static void RungeKutta2(const Real dt, const int s, Space *space, const Model *model)
{
    /*
     * Solve U(1) = LLLU = 0.0 * Un + 1.0 * LLUn
     */
    LLL(dt, 0.0, 1.0, TO, TO, TN, s, space, model);
    BoundaryConditionsAndTreatments(TN, space, model);
    /*
     * Solve U(2) = LLLU = 1.0/2.0 * Un + 1.0/2.0 * LLU(1)
     */
    LLL(dt, 1.0/2.0, 1.0/2.0, TO, TN, TO, s, space, model);
    BoundaryConditionsAndTreatments(TO, space, model);
    return;
}
static void RungeKutta3(const Real dt, const int s, Space *space, const Model *model)
{
    /*
     * Solve U(1) = LLLU = 0.0 * Un + 1.0 * LLUn
     */
    LLL(dt, 0.0, 1.0, TO, TO, TN, s, space, model);
    BoundaryConditionsAndTreatments(TN, space, model);
    /*
     * Solve U(2) = LLLU = 3.0/4.0 * Un + 1.0/4.0 * LLU(1)
     */
    LLL(dt, 3.0/4.0, 1.0/4.0, TO, TN, TM, s, space, model);
    BoundaryConditionsAndTreatments(TM, space, model);
    /*
     * Solve U(n+1) = LLLU = 1.0/3.0 * Un + 2.0/3.0 * LLU(2)
     */
    LLL(dt, 1.0/3.0, 2.0/3.0, TO, TM, TO, s, space, model);
    BoundaryConditionsAndTreatments(TO, space, model);
    return;
}
/*
 * Spatial operator computation.
 * LLLU = coeA * Un + coeB * LLU; LLU = (I + dt*L)U; LL = {LLs}; s = X, Y, Z.
 * Strategy for general coding: use s as spatial identifier, use general
 * algorithms and function pointers to unify the function and code for each
 * value of s, that is, for each spatial dimension. If a function is too
 * difficult to do a general code, then code functions for each spatial 
 * dimension individually.
 */
static void LLL(const Real dt, const Real coeA, const Real coeB, const int to, 
        const int tn, const int tm, const int p, Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    int i = 0, j = 0, k = 0; /* index with normal order */
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    Real LU[5][DIMU] = {{0.0}}; /* spatial operator */
    Real *restrict FhatR = LU[0]; /* reconstructed numerical convective flux vector */
    Real *restrict FhatL = LU[1]; /* reconstructed numerical convective flux vector */
    Real *restrict FvhatR = LU[2]; /* reconstructed numerical diffusive flux vector */
    Real *restrict FvhatL = LU[3]; /* reconstructed numerical diffusive flux vector */
    Real *restrict Phi = LU[4]; /* source vector */
    Real *temp = NULL;
    const IntVec partn = {part->n[X], part->n[Y], part->n[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const RealVec r = {dt * dd[X], dt * dd[Y], dt * dd[Z]};
    const Real rPhi = (DIMS == p) ? dt : (1.0 / 3.0) * dt;
    const int sN = (DIMS == p) ? p : p + 1;
    for (int s = (DIMS == p) ? 0 : p; s < sN; ++s) {
        for (int ks = part->np[s][Z][MIN]; ks < part->np[s][Z][MAX]; ++ks) {
            for (int js = part->np[s][Y][MIN]; js < part->np[s][Y][MAX]; ++js) {
                for (int is = part->np[s][X][MIN], state = 0; is < part->np[s][X][MAX]; ++is) {
                    switch (s) {
                        case X:
                            i = is; j = js; k = ks;
                            break;
                        case Y:
                            i = js; j = is; k = ks;
                            break;
                        case Z:
                            i = js; j = ks; k = is;
                            break;
                        default:
                            break;
                    }
                    idx = IndexNode(k, j, i, partn[Y], partn[X]);
                    if (0 != node[idx].gid) {
                        state = 0;
                        continue;
                    }
                    if (1 == state) {
                        temp = FhatL;
                        FhatL = FhatR;
                        FhatR = temp;
                        temp = FvhatL;
                        FvhatL = FvhatR;
                        FvhatR = temp;
                    } else {
                        NumericalConvectiveFlux(tn, s, k - h[s][Z], j - h[s][Y], i - h[s][X], partn, node, model, FhatL);
                        NumericalDiffusiveFlux(tn, s, k - h[s][Z], j - h[s][Y], i - h[s][X], partn, dd, node, model, FvhatL);
                        state = 1;
                    }
                    NumericalConvectiveFlux(tn, s, k, j, i, partn, node, model, FhatR);
                    NumericalDiffusiveFlux(tn, s, k, j, i, partn, dd, node, model, FvhatR);
                    SourceVector(tn, k, j, i, partn, node, model, Phi);
                    SolveOperator[model->multidim](s, coeA, coeB, node[idx].U[to], node[idx].U[tn], node[idx].U[tm], 
                            r, FhatR, FhatL, FvhatR, FvhatL, rPhi, Phi);
                }
            }
        }
    }
    return;
}
/*
 * Note: Uo, Un, and Um are all restricted pointers.
 * Under the condition that Un and Um NEVER alias each other,
 * Uo and Un may alias safely since they only read elements
 * and never modify any elements. Uo and Um may alias safely
 * since Uo only fetch the single element that Um modifies later.
 */
static void DimensionalSplitting(const int s, const Real coeA, const Real coeB, 
        const Real Uo[restrict], const Real Un[restrict], Real Um[restrict], const Real r[restrict], 
        const Real FhatR[restrict], const Real FhatL[restrict], const Real FvhatR[restrict],
        const Real FvhatL[restrict], const Real rPhi, const Real Phi[restrict])
{
    for (int dim = 0; dim < DIMU; ++dim) {
        Um[dim] = coeA * Uo[dim] + coeB * (Un[dim] - r[s] * (FhatR[dim] - FhatL[dim]) + r[s] * (FvhatR[dim] - FvhatL[dim]) + rPhi * Phi[dim]);
    }
    return;
}
static void DimensionByDimension(const int s, const Real coeA, const Real coeB, 
        const Real Uo[restrict], const Real Un[restrict], Real Um[restrict], const Real r[restrict], 
        const Real FhatR[restrict], const Real FhatL[restrict], const Real FvhatR[restrict],
        const Real FvhatL[restrict], const Real rPhi, const Real Phi[restrict])
{
    if (X == s) { /* initialize */
        for (int dim = 0; dim < DIMU; ++dim) {
            Um[dim] = coeA * Uo[dim] + coeB * (Un[dim] + rPhi * Phi[dim]);
        }
    }
    /* accumulate Ls */
    for (int dim = 0; dim < DIMU; ++dim) {
        Um[dim] = Um[dim] + coeB * (- r[s] * (FhatR[dim] - FhatL[dim]) + r[s] * (FvhatR[dim] - FvhatL[dim]));
    }
    return;
}
static void NumericalConvectiveFlux(const int tn, const int s, const int k, const int j, const int i, 
        const int partn[restrict], const Node *const node, const Model *model, Real Fhat[restrict])
{
    ReconstructConvectiveFlux[model->sScheme](tn, s, k, j, i, partn, node, model, Fhat);
    return;
}
static void NumericalDiffusiveFlux(const int tn, const int s, const int k, const int j, const int i,
        const int partn[restrict], const Real dd[restrict], const Node *const node, 
        const Model *model, Real Fvhat[restrict])
{
    const Real zero = 0.0;
    if (zero >= model->refMu) {
        memset(Fvhat, 0, DIMU * sizeof(*Fvhat));
        return;
    }
    ReconstructDiffusiveFlux[s](tn, k, j, i, partn, dd, node, model, Fvhat);
    return;
}
static void NumericalDiffusiveFluxX(const int tn, const int k, const int j, const int i,
        const int partn[restrict], const Real dd[restrict], const Node *const node, 
        const Model *model, Real Fvhat[restrict])
{
    const int idx = IndexNode(k, j, i, partn[Y], partn[X]);
    const int idxS = IndexNode(k, j - 1, i, partn[Y], partn[X]);
    const int idxN = IndexNode(k, j + 1, i, partn[Y], partn[X]);
    const int idxF = IndexNode(k - 1, j, i, partn[Y], partn[X]);
    const int idxB = IndexNode(k + 1, j, i, partn[Y], partn[X]);

    const int idxE = IndexNode(k, j, i + 1, partn[Y], partn[X]);
    const int idxSE = IndexNode(k, j - 1, i + 1, partn[Y], partn[X]);
    const int idxNE = IndexNode(k, j + 1, i + 1, partn[Y], partn[X]);
    const int idxFE = IndexNode(k - 1, j, i + 1, partn[Y], partn[X]);
    const int idxBE = IndexNode(k + 1, j, i + 1, partn[Y], partn[X]);

    const Real *restrict U = node[idx].U[tn];
    const Real u = U[1] / U[0];
    const Real v = U[2] / U[0];
    const Real w = U[3] / U[0];
    const Real T = ComputeTemperature(model->cv, U);

    U = node[idxS].U[tn];
    const Real uS = U[1] / U[0];
    const Real vS = U[2] / U[0];

    U = node[idxN].U[tn];
    const Real uN = U[1] / U[0];
    const Real vN = U[2] / U[0];

    U = node[idxF].U[tn];
    const Real uF = U[1] / U[0];
    const Real wF = U[3] / U[0];

    U = node[idxB].U[tn];
    const Real uB = U[1] / U[0];
    const Real wB = U[3] / U[0];

    U = node[idxE].U[tn];
    const Real uE = U[1] / U[0];
    const Real vE = U[2] / U[0];
    const Real wE = U[3] / U[0];
    const Real TE = ComputeTemperature(model->cv, U);

    U = node[idxSE].U[tn];
    const Real uSE = U[1] / U[0];
    const Real vSE = U[2] / U[0];

    U = node[idxNE].U[tn];
    const Real uNE = U[1] / U[0];
    const Real vNE = U[2] / U[0];

    U = node[idxFE].U[tn];
    const Real uFE = U[1] / U[0];
    const Real wFE = U[3] / U[0];

    U = node[idxBE].U[tn];
    const Real uBE = U[1] / U[0];
    const Real wBE = U[3] / U[0];

    const Real du_dx = (uE - u) * dd[X];
    const Real dv_dy = 0.25 * (vN + vNE - vS - vSE) * dd[Y];
    const Real dw_dz = 0.25 * (wB + wBE - wF - wFE) * dd[Z];
    const Real du_dy = 0.25 * (uN + uNE - uS - uSE) * dd[Y];
    const Real dv_dx = (vE - v) * dd[X];
    const Real du_dz = 0.25 * (uB + uBE - uF - uFE) * dd[Z];
    const Real dw_dx = (wE - w) * dd[X];
    const Real dT_dx = (TE - T) * dd[X];

    /* Calculate interfacial values */
    const Real uhat = 0.5 * (u + uE);
    const Real vhat = 0.5 * (v + vE);
    const Real what = 0.5 * (w + wE);
    const Real That = 0.5 * (T + TE);
    const Real mu = model->refMu * Viscosity(That * model->refT);
    const Real heatK = model->gamma * model->cv * mu / PrandtlNumber();
    const Real divV = du_dx + dv_dy + dw_dz;

    Fvhat[0] = 0.0;
    Fvhat[1] = mu * (du_dx + du_dx - (2.0/3.0) * divV);
    Fvhat[2] = mu * (du_dy + dv_dx);
    Fvhat[3] = mu * (du_dz + dw_dx);
    Fvhat[4] = heatK * dT_dx + Fvhat[1] * uhat + Fvhat[2] * vhat + Fvhat[3] * what;
    return;
}
static void NumericalDiffusiveFluxY(const int tn, const int k, const int j, const int i,
        const int partn[restrict], const Real dd[restrict], const Node *const node, 
        const Model *model, Real Fvhat[restrict])
{
    const int idx = IndexNode(k, j, i, partn[Y], partn[X]);
    const int idxW = IndexNode(k, j, i - 1, partn[Y], partn[X]);
    const int idxE = IndexNode(k, j, i + 1, partn[Y], partn[X]);
    const int idxF = IndexNode(k - 1, j, i, partn[Y], partn[X]);
    const int idxB = IndexNode(k + 1, j, i, partn[Y], partn[X]);

    const int idxN = IndexNode(k, j + 1, i, partn[Y], partn[X]);
    const int idxWN = IndexNode(k, j + 1, i - 1, partn[Y], partn[X]);
    const int idxEN = IndexNode(k, j + 1, i + 1, partn[Y], partn[X]);
    const int idxFN = IndexNode(k - 1, j + 1, i, partn[Y], partn[X]);
    const int idxBN = IndexNode(k + 1, j + 1, i, partn[Y], partn[X]);

    const Real *restrict U = node[idx].U[tn];
    const Real u = U[1] / U[0];
    const Real v = U[2] / U[0];
    const Real w = U[3] / U[0];
    const Real T = ComputeTemperature(model->cv, U);

    U = node[idxW].U[tn];
    const Real uW = U[1] / U[0];
    const Real vW = U[2] / U[0];

    U = node[idxE].U[tn];
    const Real uE = U[1] / U[0];
    const Real vE = U[2] / U[0];

    U = node[idxF].U[tn];
    const Real vF = U[2] / U[0];
    const Real wF = U[3] / U[0];

    U = node[idxB].U[tn];
    const Real vB = U[2] / U[0];
    const Real wB = U[3] / U[0];

    U = node[idxN].U[tn];
    const Real uN = U[1] / U[0];
    const Real vN = U[2] / U[0];
    const Real wN = U[3] / U[0];
    const Real TN = ComputeTemperature(model->cv, U);

    U = node[idxWN].U[tn];
    const Real uWN = U[1] / U[0];
    const Real vWN = U[2] / U[0];

    U = node[idxEN].U[tn];
    const Real uEN = U[1] / U[0];
    const Real vEN = U[2] / U[0];

    U = node[idxFN].U[tn];
    const Real vFN = U[2] / U[0];
    const Real wFN = U[3] / U[0];

    U = node[idxBN].U[tn];
    const Real vBN = U[2] / U[0];
    const Real wBN = U[3] / U[0];

    const Real dv_dx = 0.25 * (vE + vEN - vW - vWN) * dd[X];
    const Real du_dy = (uN - u) * dd[Y];
    const Real dv_dy = (vN - v) * dd[Y];
    const Real du_dx = 0.25 * (uE + uEN - uW - uWN) * dd[X];
    const Real dw_dz = 0.25 * (wB + wBN - wF - wFN) * dd[Z];
    const Real dv_dz = 0.25 * (vB + vBN - vF - vFN) * dd[Z];
    const Real dw_dy = (wN - w) * dd[Y];
    const Real dT_dy = (TN - T) * dd[Y];

    /* Calculate interfacial values */
    const Real uhat = 0.5 * (u + uN);
    const Real vhat = 0.5 * (v + vN);
    const Real what = 0.5 * (w + wN);
    const Real That = 0.5 * (T + TN);
    const Real mu = model->refMu * Viscosity(That * model->refT);
    const Real heatK = model->gamma * model->cv * mu / PrandtlNumber();
    const Real divV = du_dx + dv_dy + dw_dz;

    Fvhat[0] = 0.0;
    Fvhat[1] = mu * (dv_dx + du_dy);
    Fvhat[2] = mu * (dv_dy + dv_dy - (2.0/3.0) * divV);
    Fvhat[3] = mu * (dv_dz + dw_dy);
    Fvhat[4] = heatK * dT_dy + Fvhat[1] * uhat + Fvhat[2] * vhat + Fvhat[3] * what;
    return ;
}
static void NumericalDiffusiveFluxZ(const int tn, const int k, const int j, const int i,
        const int partn[restrict], const Real dd[restrict], const Node *const node, 
        const Model *model, Real Fvhat[restrict])
{
    const int idx = IndexNode(k, j, i, partn[Y], partn[X]);
    const int idxW = IndexNode(k, j, i - 1, partn[Y], partn[X]);
    const int idxE = IndexNode(k, j, i + 1, partn[Y], partn[X]);
    const int idxS = IndexNode(k, j - 1, i, partn[Y], partn[X]);
    const int idxN = IndexNode(k, j + 1, i, partn[Y], partn[X]);

    const int idxB = IndexNode(k + 1, j, i, partn[Y], partn[X]);
    const int idxWB = IndexNode(k + 1, j, i - 1, partn[Y], partn[X]);
    const int idxEB = IndexNode(k + 1, j, i + 1, partn[Y], partn[X]);
    const int idxSB = IndexNode(k + 1, j - 1, i, partn[Y], partn[X]);
    const int idxNB = IndexNode(k + 1, j + 1, i, partn[Y], partn[X]);

    const Real *restrict U = node[idx].U[tn];
    const Real u = U[1] / U[0];
    const Real v = U[2] / U[0];
    const Real w = U[3] / U[0];
    const Real T = ComputeTemperature(model->cv, U);

    U = node[idxW].U[tn];
    const Real uW = U[1] / U[0];
    const Real wW = U[3] / U[0];

    U = node[idxE].U[tn];
    const Real uE = U[1] / U[0];
    const Real wE = U[3] / U[0];

    U = node[idxS].U[tn];
    const Real vS = U[2] / U[0];
    const Real wS = U[3] / U[0];

    U = node[idxN].U[tn];
    const Real vN = U[2] / U[0];
    const Real wN = U[3] / U[0];

    U = node[idxB].U[tn];
    const Real uB = U[1] / U[0];
    const Real vB = U[2] / U[0];
    const Real wB = U[3] / U[0];
    const Real TB = ComputeTemperature(model->cv, U);

    U = node[idxWB].U[tn];
    const Real uWB = U[1] / U[0];
    const Real wWB = U[3] / U[0];

    U = node[idxEB].U[tn];
    const Real uEB = U[1] / U[0];
    const Real wEB = U[3] / U[0];

    U = node[idxSB].U[tn];
    const Real vSB = U[2] / U[0];
    const Real wSB = U[3] / U[0];

    U = node[idxNB].U[tn];
    const Real vNB = U[2] / U[0];
    const Real wNB = U[3] / U[0];

    const Real dw_dx = 0.25 * (wE + wEB - wW - wWB) * dd[X];
    const Real du_dz = (uB - u) * dd[Z];
    const Real dw_dy = 0.25 * (wN + wNB - wS - wSB) * dd[Y];
    const Real dv_dz = (vB - v) * dd[Z];
    const Real du_dx = 0.25 * (uE + uEB - uW - uWB) * dd[X];
    const Real dv_dy = 0.25 * (vN + vNB - vS - vSB) * dd[Y];
    const Real dw_dz = (wB - w) * dd[Z];
    const Real dT_dz = (TB - T) * dd[Z];

    /* Calculate interfacial values */
    const Real uhat = 0.5 * (u + uB);
    const Real vhat = 0.5 * (v + vB);
    const Real what = 0.5 * (w + wB);
    const Real That = 0.5 * (T + TB);
    const Real mu = model->refMu * Viscosity(That * model->refT);
    const Real heatK = model->gamma * model->cv * mu / PrandtlNumber();
    const Real divV = du_dx + dv_dy + dw_dz;

    Fvhat[0] = 0.0;
    Fvhat[1] = mu * (dw_dx + du_dz);
    Fvhat[2] = mu * (dw_dy + dv_dz);
    Fvhat[3] = mu * (dw_dz + dw_dz - (2.0/3.0) * divV);
    Fvhat[4] = heatK * dT_dz + Fvhat[1] * uhat + Fvhat[2] * vhat + Fvhat[3] * what;
    return;
}
/*
 * Source vector is splitted into three identical entities to ensure
 * consistency with spatial splitting.
 */
static void SourceVector(const int tn, const int k, const int j, const int i,
        const int partn[restrict], const Node *const node, const Model *model, Real Phi[restrict])
{
    if (0 == model->sState) {
        memset(Phi, 0, DIMU * sizeof(*Phi));
        return;
    }
    const int idx = IndexNode(k, j, i, partn[Y], partn[X]);
    const Real *restrict U = node[idx].U[tn];
    const RealVec V = {U[1] / U[0], U[2] / U[0], U[3] / U[0]};
    const RealVec fb = {U[0] * model->g[X], U[0] * model->g[Y], U[0] * model->g[Z]};
    Phi[0] = 0.0;
    Phi[1] = fb[X];
    Phi[2] = fb[Y];
    Phi[3] = fb[Z];
    Phi[4] = Dot(fb, V);
    return;
}
/* a good practice: end file with a newline */

