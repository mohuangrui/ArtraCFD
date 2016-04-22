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
typedef int (*NumericalFluxReconstructor)(const int, Real [], const Real, const int,
        const int, const int, const Real *, const Space *, const Model *);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int DifferentialOperatorSplitting(Field *, const Space *, const Model *,
        const Partition *, const Geometry *, const Real);
static int RungeKutta(const int, Field *, const Space *, const Model *,
        const Partition *, const Geometry *, const Real);
static int LL(const int, Real *, const Real *, const Space *, const Model *,
        const Partition *, const Real);
static int NumericalConvectiveFlux(const int, Real [], const Real, const int,
        const int, const int, const Real *, const Space *, const Model *);
static int DiffusiveFluxGradient(const int, Real [], const int, const int,
        const int, const Real *, const Space *, const Model *);
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
int FluidDynamics(Field *field, const Space *space, const Model *model, 
        const Partition *part, const Geometry *geometry, const Real dt)
{
    DifferentialOperatorSplitting(field, space, model, part, geometry, dt);
    return 0;
}
static int DifferentialOperatorSplitting(Field *field, const Space *space, const Model *model,
        const Partition *part, const Geometry *geometry, const Real dt)
{
    RungeKutta(Z, field, space, model, part, geometry, 0.5 * dt);
    RungeKutta(Y, field, space, model, part, geometry, 0.5 * dt);
    RungeKutta(X, field, space, model, part, geometry, 0.5 * dt);
    RungeKutta(X, field, space, model, part, geometry, 0.5 * dt);
    RungeKutta(Y, field, space, model, part, geometry, 0.5 * dt);
    RungeKutta(Z, field, space, model, part, geometry, 0.5 * dt);
    return 0;
}
static int RungeKutta(const int s, Field *field, const Space *space, const Model *model,
        const Partition *part, const Geometry *geometry, const Real dt)
{
    Real *exchanger = field->U;
    /*
     * First, save the full value of current field, since this value is required
     * for the calculation of field data at intermediate stage as well as n+1.
     * Operation can be achieved by a single loop since all data are stored by
     * linear arrays.
     */
    for (int idx = 0; idx < (space->nMax * DIMU); ++idx) {
        field->Un[idx] = field->U[idx];
    }
    /*
     * Then solve LLU = (I + dt*L)U.
     */
    /*
     * When exchange a large bunch of data between two storage space, such as
     * arrays, if there is no new data generation but just data exchange and 
     * update, then the rational way is to exchange the address value that
     * their pointer point to rather than values of data entries.
     */
    LL(s, field->Uswap, field->U, space, model, part, dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part, geometry);
    exchanger = field->U; /* preserve the address of U */
    field->U = field->Uswap; /* update flow field */
    field->Uswap = exchanger; /* regain the used space as new space */
    /*
     * Now solve the LLU = (I + dt*L)U based on updated field for another time step.
     */
    LL(s, field->Uswap, field->U, space, model, part, dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part, geometry);
    exchanger = field->U; /* preserve the address of U */
    field->U = field->Uswap; /* update flow field */
    field->Uswap = exchanger; /* regain the used space as new space */
    /*
     * Calculate the intermediate stage based on the newest updated data and
     * the original field data which is stored at first. No new storage space
     * need to be introduced since all calculations are based on a single space
     * and do not require neighbours.
     */
    const Real coeA = 3.0 / 4.0; /* caution! float operation required! */
    const Real coeB = 1.0 / 4.0; /* caution! float operation required! */
    for (int idx = 0; idx < (space->nMax * DIMU); ++idx) {
        field->U[idx] = coeA * field->Un[idx] + coeB * field->U[idx];
    }
    /*
     * Now solve LLU = (I + dt*L)U based on the intermediate field.
     */
    LL(s, field->Uswap, field->U, space, model, part, dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part, geometry);
    exchanger = field->U; /* preserve the address of U */
    field->U = field->Uswap; /* update flow field */
    field->Uswap = exchanger; /* regain the used space as new space */
    /*
     * Calculate field data at n+1
     */
    const Real coeAA = 1.0 / 3.0; /* caution! float operation required! */
    const Real coeBB = 2.0 / 3.0; /* caution! float operation required! */
    for (int idx = 0; idx < (space->nMax * DIMU); ++idx) {
        field->U[idx] = coeAA * field->Un[idx] + coeBB * field->U[idx];
    }
    return 0;
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
static int LL(const int s, Real *U, const Real *Un, const Space *space,
        const Model *model, const Partition *part, const Real dt)
{
    Real Fhat[DIMU] = {0.0}; /* reconstructed numerical convective flux vector */
    Real Fhath[DIMU] = {0.0}; /* reconstructed numerical convective flux vector at neighbour */
    Real gradG[DIMU] = {0.0}; /* spatial gradient of diffusive flux vector */
    const int h[DIMS][DIMS] = {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}; /* direction indicator */
    const Real r[DIMS] = {dt * space->ddx, dt * space->ddy, dt * space->ddz};
    int idx = 0; /* linear array index math variable */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (FLUID != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                NumericalConvectiveFlux(s, Fhat, r[s], k, j, i, Un, space, model);
                NumericalConvectiveFlux(s, Fhath, r[s], k + h[s][Z], j + h[s][Y], i + h[s][X], Un, space, model);
                DiffusiveFluxGradient(s, gradG, k, j, i, Un, space, model);
                idx = idx * DIMU; /* change idx to field variable */
                for (int dim = 0; dim < DIMU; ++dim) {
                    /* conservative discretization for convective flux */
                    U[idx+dim] = Un[idx+dim] - r[s] * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int NumericalConvectiveFlux(const int s, Real Fhat[], const Real r,
        const int k, const int j, const int i,
        const Real *U, const Space *space, const Model *model)
{
    NumericalFluxReconstructor ReconstructNumericalFlux[2] = {
        TVD,
        WENO
    };
    ReconstructNumericalFlux[model->scheme](s, Fhat, r, k, j, i, U, space, model);
    return 0;
}
/*
 * Generally the diffusive terms will only be discretized by central
 * difference scheme. However, for outermost layer of interior nodes, the
 * central scheme of them require the diffusive flux vector at boundaries, 
 * because the corners of current computational domain haven't set with any
 * values, the diffusive flux vector at boundaries can not be interpolated
 * with central scheme, therefore, need to change to forward or backward. 
 * This situation is the same to interior ghost nodes the central
 * difference scheme can't be applied because of lacking stencil. Thus,
 * they also need to be identified and modified.
 */
static int DiffusiveFluxGradient(const int s, Real gradG[],
        const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real Gl[DIMU] = {0.0}; /* diffusive flux vector at left */
    Real Gr[DIMU] = {0.0}; /* diffusive flux vector at right */
    Real dL[DIMS] = {space->ddx, space->ddy, space->ddz}; /* reciprocal of differencing distance */
    Real dCoe = 0;
    /* default is central scheme */
    int hl[DIMS][DIMS] = {{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}; /* left offset */
    int hr[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* right offset */
    if (0 == (k + hl[s][Z] - space->ng) * (j + hl[s][Y] - space->ng) * 
            (i + hl[s][X] - space->ng)) { /* if left at boundary, no offset */
        hl[s][Z] = 0;
        hl[s][Y] = 0;
        hl[s][X] = 0;
    }
    if (0 == (k + hr[s][Z] - space->ng - space->nz + 1) * (j + hr[s][Y] - space->ng - space->ny + 1) * 
            (i + hr[s][X] - space->ng - space->nx + 1)) { /* if right at boundary, no offset */
        hr[s][Z] = 0;
        hr[s][Y] = 0;
        hr[s][X] = 0;
    }
    /* check ghost */
    const int idxl = IndexMath(k + hl[s][Z], j + hl[s][Y], i + hl[s][X], space);
    const int idxr = IndexMath(k + hr[s][Z], j + hr[s][Y], i + hr[s][X], space);
    if (OFFSET <= space->nodeFlag[idxl]) { /* if left is ghost, no offset */
        hl[s][Z] = 0;
        hl[s][Y] = 0;
        hl[s][X] = 0;
    }
    if (OFFSET <= space->nodeFlag[idxr]) { /* if right is ghost, no offset */
        hr[s][Z] = 0;
        hr[s][Y] = 0;
        hr[s][X] = 0;
    }
    dCoe = hr[s][Z] - hl[s][Z] + hr[s][Y] - hl[s][Y] + hr[s][X] - hr[s][X];
    if (0 != dCoe) { /* only do calculation when offset exists */
        dL[s] = dL[s] / dCoe;
        DiffusiveFlux(s, Gl, k + hl[s][Z], j + hl[s][Y], i + hl[s][X], U, space, model);
        DiffusiveFlux(s, Gr, k + hr[s][Z], j + hr[s][Y], i + hr[s][X], U, space, model);
    }
    for (int row = 0; row < DIMU; ++row) {
        gradG[row] = dL[s] * (Gr[row] - Gl[row]);
    }
    return 0;
}
/* a good practice: end file with a newline */

