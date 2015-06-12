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
#include "fluid_dynamics.h"
#include <stdio.h> /* standard library for input and output */
#include "tvd.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int RungeKutta(Field *, const Space *, const Model *, const Partition *,
        Geometry *, const Real);
static int AlgebraicOperatorSplitting(Field *, const Space *, const Model *, 
        const Partition *, const Geometry *, const Real);
static int Lz(Real *U, const Real *Un, const Space *space, 
        const Partition *part, const Flow *flow, const Real dt);
static int Ly(Real *U, const Real *Un, const Space *space, 
        const Partition *part, const Flow *flow, const Real dt);
static int Lx(Real *U, const Real *Un, const Space *space, 
        const Partition *part, const Flow *flow, const Real dt);
static int ComputeReconstructedFluxZ(
        Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeReconstructedFluxY(
        Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
static int ComputeReconstructedFluxX(
        Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Flow *flow, const Real dt);
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
 *       - spactial discretization (convective fluxes + diffusive fluxes)
 */
int FluidDynamics(Field *field, const Space *space, const Model *model, 
        const Partition *part, Geometry *geometry, const Real dt)
{
    RungeKutta(field, space, model, part, geometry, dt);
}
static int RungeKutta(Field *field, const Space *space, const Model *model,
        const Partition *part, Geometry *geometry, const Real dt)
{
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
     * Then solve LLU = (I + dt*L)U by algebraic operator splitting, updated data 
     * will be stored in the same storage space of inputed data.
     */
    AlgebraicOperatorSplitting(field->U, field->Uswap, space, model, part, geometry, dt);
    /*
     * Now solve the LLU = (I + dt*L)U based on updated field for another time step.
     */
    AlgebraicOperatorSplitting(field->U, field->Uswap, space, model, part, geometry, dt);
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
    AlgebraicOperatorSplitting(field, space, model, part, geometry, dt);
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
static int AlgebraicOperatorSplitting(Field *field, const Space *space,
        const Model *model, const Partition *part, const Geometry *geometry, const Real dt)
{
    Real *exchanger = U;
    /*
     * When exchange a large bunch of data between two storage space, such as
     * arrays, if there is no new data generation but just data exchange and 
     * update, then the rational way is to exchange the address value that
     * their pointer point to rather than values of data entries.
     */
    LLz(Uswap, U, space, model, part, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, geometry, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Ly(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, geometry, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Lx(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, geometry, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Lx(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, geometry, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Ly(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, geometry, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */

    Lz(Uswap, U, space, part, flow, 0.5 * dt);
    BoundaryCondtionsAndTreatments(Uswap, space, geometry, part, flow);
    exchanger = U; /* preserve the address of U */
    U = Uswap; /* update flow field */
    Uswap = exchanger; /* regain the used space as new space */
    /*
     * NOTICE: the exchange operations must be even times, otherwise, the
     * storage space that U originally pointed to will not store the newest
     * updated field data after computation.
     */
    return 0;
}
static int Lz(Real *U, const Real *Un, const Space *space, const Partition *part, const Flow *flow, const Real dt)
{
    Real Fhat[DIMU] = {0.0}; /* reconstructed flux vector */
    Real Fhath[DIMU] = {0.0}; /* reconstructed flux vector at neighbour */
    Real gradG[DIMU] = {0.0}; /* spatial gradient of viscous flux vector */
    int idx = 0; /* linear array index math variable */
    const Real r = dt * space->ddz;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (FLUID != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                ComputeReconstructedFluxZ(Fhat, k, j, i, Un, space, flow, dt);
                ComputeReconstructedFluxZ(Fhath, k - 1, j, i, Un, space, flow, dt);
                ComputeViscousFluxGradientZ(gradG, k, j, i, Un, space, flow);
                idx = idx * DIMU; /* change idx to field variable */
                for (int dim = 0; dim < DIMU; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int Ly(Real *U, const Real *Un, const Space *space, const Partition *part, const Flow *flow, const Real dt)
{
    Real Fhat[DIMU] = {0.0}; /* reconstructed flux vector */
    Real Fhath[DIMU] = {0.0}; /* reconstructed flux vector at neighbour */
    Real gradG[DIMU] = {0.0}; /* spatial gradient of viscous flux vector */
    int idx = 0; /* linear array index math variable */
    const Real r = dt * space->ddy;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (FLUID != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                ComputeReconstructedFluxY(Fhat, k, j, i, Un, space, flow, dt);
                ComputeReconstructedFluxY(Fhath, k, j - 1, i, Un, space, flow, dt);
                ComputeViscousFluxGradientY(gradG, k, j, i, Un, space, flow);
                idx = idx * DIMU; /* change idx to field variable */
                for (int dim = 0; dim < DIMU; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int Lx(Real *U, const Real *Un, const Space *space, const Partition *part, const Flow *flow, const Real dt)
{
    Real Fhat[DIMU] = {0.0}; /* reconstructed flux vector */
    Real Fhath[DIMU] = {0.0}; /* reconstructed flux vector at neighbour */
    Real gradG[DIMU] = {0.0}; /* spatial gradient of viscous flux vector */
    int idx = 0; /* linear array index math variable */
    const Real r = dt * space->ddx;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (FLUID != space->nodeFlag[idx]) { /* it's not a fluid */
                    continue;
                }
                ComputeReconstructedFluxX(Fhat, k, j, i, Un, space, flow, dt);
                ComputeReconstructedFluxX(Fhath, k, j, i - 1, Un, space, flow, dt);
                ComputeViscousFluxGradientX(gradG, k, j, i, Un, space, flow);
                idx = idx * DIMU; /* change idx to field variable */
                for (int dim = 0; dim < DIMU; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
/* a good practice: end file with a newline */

