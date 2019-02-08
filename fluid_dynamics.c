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
#include "convective_flux.h"
#include "diffusive_flux.h"
#include "source_term.h"
#include "boundary_treatment.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Function Pointers
 ****************************************************************************/
typedef void (*TimeIntegrator)(const Real, const int, Space *, const Model *);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void DiscretizeTime(const Real, const int, Space *, const Model *);
static void RungeKutta2(const Real, const int, Space *, const Model *);
static void RungeKutta3(const Real, const int, Space *, const Model *);
static void LLLU(const Real, const Real, const Real, const int,
        const int, const int, const int, Space *, const Model *);
static void LU(const Real [restrict], const Real [restrict],
        const Real [restrict], const Real [restrict], Real [restrict]);
static void SolveOperator(const int, const int, const Real, const Real,
        const Real [restrict], const Real [restrict], Real [restrict], const Real,
        const Real [restrict]);
/****************************************************************************
 * Global Variables Definition with Private Scope
 ****************************************************************************/
static TimeIntegrator IntegrateTime[2] = {
    RungeKutta2,
    RungeKutta3};
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * dU/dt = LU = LxU + LyU + LzU + Phi(U)
 * Time and space discretizations are implemented under the method of lines.
 * Operator splitting is used to treat equation with source terms.
 * Multi-dimensionality is addressed by two approaches
 *   a) - operator splitting
 *   b) - operator-by-operator approximation
 */
void EvolveFluidDynamics(const Real dt, Space *space, const Model *model)
{
    if (0 != model->sState) {
        DiscretizeTime(0.5 * dt, PHI, space, model);
    }
    switch (model->multidim) {
        case OPTSPLIT:
            switch (space->part.collapse) {
                case COLLAPSEN:
                    DiscretizeTime(0.5 * dt, Z, space, model);
                    DiscretizeTime(0.5 * dt, Y, space, model);
                    DiscretizeTime(0.5 * dt, X, space, model);
                    DiscretizeTime(0.5 * dt, X, space, model);
                    DiscretizeTime(0.5 * dt, Y, space, model);
                    DiscretizeTime(0.5 * dt, Z, space, model);
                    break;
                case COLLAPSEX:
                    DiscretizeTime(0.5 * dt, Z, space, model);
                    DiscretizeTime(0.5 * dt, Y, space, model);
                    DiscretizeTime(0.5 * dt, Y, space, model);
                    DiscretizeTime(0.5 * dt, Z, space, model);
                    break;
                case COLLAPSEY:
                    DiscretizeTime(0.5 * dt, Z, space, model);
                    DiscretizeTime(0.5 * dt, X, space, model);
                    DiscretizeTime(0.5 * dt, X, space, model);
                    DiscretizeTime(0.5 * dt, Z, space, model);
                    break;
                case COLLAPSEZ:
                    DiscretizeTime(0.5 * dt, Y, space, model);
                    DiscretizeTime(0.5 * dt, X, space, model);
                    DiscretizeTime(0.5 * dt, X, space, model);
                    DiscretizeTime(0.5 * dt, Y, space, model);
                    break;
                case COLLAPSEXY:
                    DiscretizeTime(0.5 * dt, Z, space, model);
                    DiscretizeTime(0.5 * dt, Z, space, model);
                    break;
                case COLLAPSEXZ:
                    DiscretizeTime(0.5 * dt, Y, space, model);
                    DiscretizeTime(0.5 * dt, Y, space, model);
                    break;
                case COLLAPSEYZ:
                    DiscretizeTime(0.5 * dt, X, space, model);
                    DiscretizeTime(0.5 * dt, X, space, model);
                    break;
                default:
                    break;
            }
            break;
        case OPTBYOPT:
            DiscretizeTime(0.5 * dt, DIMS, space, model);
            DiscretizeTime(0.5 * dt, DIMS, space, model);
            break;
        default:
            break;
    }
    if (0 != model->sState) {
        DiscretizeTime(0.5 * dt, PHI, space, model);
    }
    return;
}
/*
 * dU/dt = LU
 * Computation must start from TO data space and end with TO data space.
 */
static void DiscretizeTime(const Real dt, const int s, Space *space, const Model *model)
{
    IntegrateTime[model->tScheme](dt, s, space, model);
    return;
}
static void RungeKutta2(const Real dt, const int s, Space *space, const Model *model)
{
    /* solve U1 = LLLU = 0.0 * Un + 1.0 * LLUn */
    LLLU(dt, 0.0, 1.0, TO, TO, TN, s, space, model);
    TreatBoundary(TN, space, model);
    /* solve U(n+1) = LLLU = 1.0/2.0 * Un + 1.0/2.0 * LLU1 */
    LLLU(dt, 1.0/2.0, 1.0/2.0, TO, TN, TO, s, space, model);
    TreatBoundary(TO, space, model);
    return;
}
static void RungeKutta3(const Real dt, const int s, Space *space, const Model *model)
{
    /* solve U1 = LLLU = 0.0 * Un + 1.0 * LLUn */
    LLLU(dt, 0.0, 1.0, TO, TO, TN, s, space, model);
    TreatBoundary(TN, space, model);
    /* solve U2 = LLLU = 3.0/4.0 * Un + 1.0/4.0 * LLU1 */
    LLLU(dt, 3.0/4.0, 1.0/4.0, TO, TN, TM, s, space, model);
    TreatBoundary(TM, space, model);
    /* solve U(n+1) = LLLU = 1.0/3.0 * Un + 2.0/3.0 * LLU2 */
    LLLU(dt, 1.0/3.0, 2.0/3.0, TO, TM, TO, s, space, model);
    TreatBoundary(TO, space, model);
    return;
}
/*
 * Spatial operator computation.
 * LLLU = coeA * Un + coeB * LLU; LLU = (I + dt*L)U; L = {Ls, phi}; s = X, Y, Z.
 * Strategy for general coding: use p as operator identifier, use general
 * algorithms and function pointers to unify the function and code for each
 * value of p. If a function is too difficult to do general coding, then code
 * functions for each operator individually.
 */
static void LLLU(const Real dt, const Real coeA, const Real coeB, const int to,
        const int tn, const int tm, const int p, Space *space, const Model *model)
{
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    int i = 0, j = 0, k = 0; /* index with normal order */
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    Real RHS[5][DIMU] = {{0.0}}; /* spatial operator */
    Real *restrict FhatR = RHS[0]; /* reconstructed numerical convective flux vector */
    Real *restrict FhatL = RHS[1]; /* reconstructed numerical convective flux vector */
    Real *restrict FvhatR = RHS[2]; /* reconstructed numerical diffusive flux vector */
    Real *restrict FvhatL = RHS[3]; /* reconstructed numerical diffusive flux vector */
    Real *restrict Phi = RHS[4]; /* right hand side vector */
    Real *temp = NULL;
    const IntVec partn = {part->n[X], part->n[Y], part->n[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const RealVec r = {dt * dd[X], dt * dd[Y], dt * dd[Z]};
    int s = 0, sN = 0; /* space sweep control for the operator p */
    switch (p) {
        case PHI: /* source term */
            s = 0; sN = s + 1;
            break;
        case DIMS: /* all spatial operators */
            s = 0; sN = DIMS;
            break;
        default: /* individual spatial operator */
            s = p; sN = s + 1;
            break;
    }
    /* space sweep with dimension priority */
    for (; s < sN; ++s) {
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
                    if (0 != node[idx].did) {
                        state = 0; /* mark domain change and boundary occurrence */
                        continue;
                    }
                    switch (p) {
                        case PHI:
                            ComputePhi(tn, k, j, i, partn, node, model, Phi);
                            SolveOperator(OPTSPLIT, s, coeA, coeB, node[idx].U[to], node[idx].U[tn], node[idx].U[tm], dt, Phi);
                            continue;
                        default:
                            break;
                    }
                    switch (state) {
                        case 1: /* inherit numerical flux from the previous node */
                            temp = FhatL;
                            FhatL = FhatR;
                            FhatR = temp;
                            temp = FvhatL;
                            FvhatL = FvhatR;
                            FvhatR = temp;
                            break;
                        default: /* compute numerical flux at left interface */
                            ComputeFhat(tn, s, k - h[s][Z], j - h[s][Y], i - h[s][X], partn, node, model, FhatL);
                            ComputeFvhat(tn, s, k - h[s][Z], j - h[s][Y], i - h[s][X], partn, dd, node, model, FvhatL);
                            state = 1;
                            break;
                    }
                    ComputeFhat(tn, s, k, j, i, partn, node, model, FhatR);
                    ComputeFvhat(tn, s, k, j, i, partn, dd, node, model, FvhatR);
                    LU(FhatR, FhatL, FvhatR, FvhatL, Phi);
                    SolveOperator(model->multidim, s, coeA, coeB, node[idx].U[to], node[idx].U[tn], node[idx].U[tm], r[s], Phi);
                }
            }
        }
    }
    return;
}
static void LU(const Real FhatR[restrict], const Real FhatL[restrict],
        const Real FvhatR[restrict], const Real FvhatL[restrict], Real Phi[restrict])
{
    for (int n = 0; n < DIMU; ++n) {
        Phi[n] = FhatL[n] - FhatR[n] + FvhatR[n] - FvhatL[n];
    }
    return;
}
/*
 * Solve the solution operator for time integration.
 * Note: Uo, Un, and Um are all restricted pointers. Under the condition that
 * Un and Um NEVER alias each other, Uo and Un may alias safely since they only
 * read elements and never modify any elements. Uo and Um may alias safely
 * since Uo only fetch the single element that Um modifies later.
 */
static void SolveOperator(const int p, const int s, const Real coeA, const Real coeB,
        const Real Uo[restrict], const Real Un[restrict], Real Um[restrict], const Real r,
        const Real Phi[restrict])
{
    /* accumulation step for operator-by-operator approximation */
    if ((OPTBYOPT == p) && (X != s)) {
        for (int n = 0; n < DIMU; ++n) {
            Um[n] = Um[n] + coeB * r * Phi[n];
        }
        return;
    }
    /* solve step for the solution operator */
    for (int n = 0; n < DIMU; ++n) {
        Um[n] = coeA * Uo[n] + coeB * (Un[n] + r * Phi[n]);
    }
    return;
}
/* a good practice: end file with a newline */

