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
#include "weno.h"
#include "boundary_treatment.h"
#include "cfd_commons.h"
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
typedef void (*ConvectiveFluxReconstructor)(const int, const int, const int,
        const int, const int, const int [restrict], const Node *, 
        const Model *, Real [restrict]);
typedef void (*DiffusiveFluxReconstructor)(const int, const int, const int, 
        const int, const int [restrict], const Real [restrict], const Node *, 
        const Model *, Real [restrict]);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void LL(const Real, const int, const int, const int, Space *, const Model *);
static void NumericalConvectiveFlux(const int, const int, const int, const int,
        const int, const int [restrict], const Node *, const Model *, Real [restrict]);
static void NumericalDiffusiveFlux(const int, const int, const int, const int, 
        const int, const int [restrict], const Real [restrict], const Node *, 
        const Model *, Real [restrict]);
static void NumericalDiffusiveFluxX(const int, const int, const int, 
        const int, const int [restrict], const Real [restrict], const Node *, 
        const Model *, Real [restrict]);
static void NumericalDiffusiveFluxY(const int, const int, const int, 
        const int, const int [restrict], const Real [restrict], const Node *, 
        const Model *, Real [restrict]);
static void NumericalDiffusiveFluxZ(const int, const int, const int, 
        const int, const int [restrict], const Real [restrict], const Node *, 
        const Model *, Real [restrict]);
static Real Viscosity(const Real);
static Real PrandtlNumber(void);
/****************************************************************************
 * Global Variables Definition with Private Scope
 ****************************************************************************/
static ConvectiveFluxReconstructor ReconstructConvectiveFlux[1] = {
    WENO};
static DiffusiveFluxReconstructor ReconstructDiffusiveFlux[DIMS] = {
    NumericalDiffusiveFluxX,
    NumericalDiffusiveFluxY,
    NumericalDiffusiveFluxZ};
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Two approaches for the numerical solution of governing equations:
 *   a)  - differential operator splitting
 *       - temporal discretization
 *       - spatial discretization (convective fluxes + diffusive fluxes)
 *   b)  - temporal discretization
 *       - algebraic operator splitting
 *       - spatial discretization (convective fluxes + diffusive fluxes)
 *   a) shows more accurate results according to tests on Riemann problems.
 */
int FluidDynamics(Space *space, const Model *model, const Real dt)
{
    DifferentialOperatorSplitting(space, model, dt);
    return 0;
}
static int DifferentialOperatorSplitting(Space *space, const Model *model, const Real dt)
{
    RungeKutta(Z, space, model, 0.5 * dt);
    RungeKutta(Y, space, model, 0.5 * dt);
    RungeKutta(X, space, model, 0.5 * dt);
    RungeKutta(X, space, model, 0.5 * dt);
    RungeKutta(Y, space, model, 0.5 * dt);
    RungeKutta(Z, space, model, 0.5 * dt);
    return 0;
}
static void RungeKutta(const Real dt, const int s, Space *space, const Model *model)
{
    /*
     * Solve the LLU = (I + dt*L)U
     */
    LL(dt, s, C, CN, space, model);
    BoundaryCondtionsAndTreatments(CN, space, model);
    /*
     * Now solve the LLU = (I + dt*L)U based on updated field for another time step.
     */
    LL(dt, s, CN, CM, space, model);
    BoundaryCondtionsAndTreatments(CN, space, model);
    /*
     * Calculate the intermediate stage based on the newest updated data and
     * the original field data.
     */
    const Real coeA = 3.0 / 4.0; /* caution! float operation required! */
    const Real coeB = 1.0 / 4.0; /* caution! float operation required! */
    /*
     * Now solve LLU = (I + dt*L)U based on the intermediate field.
     */
    LL(s, field->Uswap, field->U, space, model, part, dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part, geometry);
    /*
     * Calculate field data at n+1
     */
    const Real coeAA = 1.0 / 3.0; /* caution! float operation required! */
    const Real coeBB = 2.0 / 3.0; /* caution! float operation required! */
    for (int idx = 0; idx < (space->nMax * DIMU); ++idx) {
        field->U[idx] = coeAA * field->Un[idx] + coeBB * field->U[idx];
    }
    return;
}
/*
 * Spatial operator computation.
 * LLU = (I + dt*L)U; LL = {LLs}; s = 0, 1, 2 for x, y, z, respectively.
 * Strategy for general coding: use s as spatial identifier, use general
 * algorithms and function pointers to unify the function and code for each
 * value of s, that is, for each spatial dimension. If a function is too
 * difficult to do a general code, then code functions for each spatial 
 * dimension individually.
 */
static void LL(const Real dt, const int s, const int tn, const int tc, Space *space,
        const Model *model)
{
    int idx = 0; /* linear array index math variable */
    Node *node = space->node;
    const Real *restrict Un = NULL;
    Real *restrict Uc = NULL;
    const Partition *restrict part = &(space->part);
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    Real FhatR[DIMU] = {0.0}; /* reconstructed numerical convective flux vector */
    Real FhatL[DIMU] = {0.0}; /* reconstructed numerical convective flux vector */
    Real FvhatR[DIMU] = {0.0}; /* reconstructed numerical diffusive flux vector */
    Real FvhatL[DIMU] = {0.0}; /* reconstructed numerical diffusive flux vector */
    const int partn[DIMS] = {part->n[X], part->n[Y], part->n[Z]};
    const Real dd[DIMS] = {part->dd[X], part->dd[Y], part->dd[Z]};
    const Real r[DIMS] = {dt * dd[X], dt * dd[Y], dt * dd[Z]};
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, partn[Y], partn[X]);
                if (FLUID != node[idx].geoID) {
                    continue;
                }
                Un = node[idx].U[tn];
                Uc = node[idx].U[tc];
                NumericalConvectiveFlux(s, tn, k, j, i, partn, node, model, FhatR);
                NumericalConvectiveFlux(s, tn, k - h[s][Z], j - h[s][Y], i - h[s][X], partn, node, model, FhatL);
                if (0 < model->refMu) {
                    NumericalDiffusiveFlux(s, tn, k, j, i, partn, dd, node, model, FvhatR);
                    NumericalDiffusiveFlux(s, tn, k - h[s][Z], j - h[s][Y], i - h[s][X], partn, dd, node, model, FvhatL);
                }
                for (int dim = 0; dim < DIMU; ++dim) {
                    /* conservative discretization for fluxes */
                    Uc[dim] = Un[dim] - r[s] * (FhatR[dim] - FhatL[dim]) + r[s] * (FvhatR[dim] - FvhatL[dim]);
                }
            }
        }
    }
    return;
}
static void NumericalConvectiveFlux(const int s, const int tn, const int k, const int j, const int i, 
        const int partn[restrict], const Node *node, const Model *model, Real Fhat[restrict])
{
    ReconstructConvectiveFlux[model->scheme](s, tn, k, j, i, partn, node, model, Fhat);
    return;
}
static void NumericalDiffusiveFlux(const int s, const int tn, const int k, const int j, 
        const int i, const int partn[restrict], const Real dd[restrict], const Node *node, 
        const Model *model, Real Fvhat[restrict])
{
    ReconstructDiffusiveFlux[s](tn, k, j, i, partn, dd, node, model, Fvhat);
    return;
}
static void NumericalDiffusiveFluxX(const int tn, const int k, const int j, 
        const int i, const int partn[restrict], const Real dd[restrict], const Node *node, 
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

    Fvhat[0] = 0;
    Fvhat[1] = mu * (2.0 * du_dx - (2.0/3.0) * divV);
    Fvhat[2] = mu * (du_dy + dv_dx);
    Fvhat[3] = mu * (du_dz + dw_dx);
    Fvhat[4] = heatK * dT_dx + Fvhat[1] * uhat + Fvhat[2] * vhat + Fvhat[3] * what;
    return;
}
static void NumericalDiffusiveFluxY(const int tn, const int k, const int j, 
        const int i, const int partn[restrict], const Real dd[restrict], const Node *node, 
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

    Fvhat[0] = 0;
    Fvhat[1] = mu * (dv_dx + du_dy);
    Fvhat[2] = mu * (2.0 * dv_dy - (2.0/3.0) * divV);
    Fvhat[3] = mu * (dv_dz + dw_dy);
    Fvhat[4] = heatK * dT_dy + Fvhat[1] * uhat + Fvhat[2] * vhat + Fvhat[3] * what;
    return ;
}
static void NumericalDiffusiveFluxZ(const int tn, const int k, const int j, 
        const int i, const int partn[restrict], const Real dd[restrict], const Node *node, 
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

    Fvhat[0] = 0;
    Fvhat[1] = mu * (dw_dx + du_dz);
    Fvhat[2] = mu * (dw_dy + dv_dz);
    Fvhat[3] = mu * (2.0 * dw_dz - (2.0/3.0) * divV);
    Fvhat[4] = heatK * dT_dz + Fvhat[1] * uhat + Fvhat[2] * vhat + Fvhat[3] * what;
    return;
}
static Real Viscosity(const Real T)
{
    return 1.458e-6 * pow(T, 1.5) / (T + 110.4); /* Sutherland's law */
}
static Real PrandtlNumber(void)
{
    return 0.71; /* air */
}
/* a good practice: end file with a newline */

