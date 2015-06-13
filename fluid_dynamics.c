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
#include <math.h> /* common mathematical functions */
#include "tvd.h"
#include "boundary_treatment.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int RungeKutta(Field *, const Space *, const Model *, const Partition *,
        const Geometry *, const Real);
static int AlgebraicOperatorSplitting(Field *, const Space *, const Model *, 
        const Partition *, const Geometry *, const Real);
static int LLz(Real *, const Real *, const Space *, const Model *,
        const Partition *, const Real);
static int LLy(Real *, const Real *, const Space *, const Model *,
        const Partition *, const Real);
static int LLx(Real *, const Real *, const Space *, const Model *,
        const Partition *, const Real);
static int ComputeReconstructedFluxZ(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *, const Real);
static int ComputeReconstructedFluxY(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *, const Real);
static int ComputeReconstructedFluxX(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *, const Real);
static int ComputeViscousFluxGradientZ(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
static int ComputeViscousFluxGradientY(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
static int ComputeViscousFluxGradientX(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
static int ComputeViscousFluxZ(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
static int ComputeViscousFluxY(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
static int ComputeViscousFluxX(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
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
        const Partition *part, const Geometry *geometry, const Real dt)
{
    RungeKutta(field, space, model, part, geometry, dt);
    return 0;
}
static int RungeKutta(Field *field, const Space *space, const Model *model,
        const Partition *part, const Geometry *geometry, const Real dt)
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
    AlgebraicOperatorSplitting(field, space, model, part, geometry, dt);
    /*
     * Now solve the LLU = (I + dt*L)U based on updated field for another time step.
     */
    AlgebraicOperatorSplitting(field, space, model, part, geometry, dt);
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
static int AlgebraicOperatorSplitting(Field *field, const Space *space, const Model *model,
        const Partition *part, const Geometry *geometry, const Real dt)
{
    Real *exchanger = field->U;
    /*
     * When exchange a large bunch of data between two storage space, such as
     * arrays, if there is no new data generation but just data exchange and 
     * update, then the rational way is to exchange the address value that
     * their pointer point to rather than values of data entries.
     */
    LLz(field->Uswap, field->U, space, model, part, 0.5 * dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part, geometry);
    exchanger = field->U; /* preserve the address of U */
    field->U = field->Uswap; /* update flow field */
    field->Uswap = exchanger; /* regain the used space as new space */

    LLy(field->Uswap, field->U, space, model, part, 0.5 * dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part,geometry);
    exchanger = field->U; /* preserve the address of U */
    field->U = field->Uswap; /* update flow field */
    field->Uswap = exchanger; /* regain the used space as new space */

    LLx(field->Uswap, field->U, space, model, part, 0.5 * dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part, geometry);
    exchanger = field->U; /* preserve the address of U */
    field->U = field->Uswap; /* update flow field */
    field->Uswap = exchanger; /* regain the used space as new space */

    LLx(field->Uswap, field->U, space, model, part, 0.5 * dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part, geometry);
    exchanger = field->U; /* preserve the address of U */
    field->U = field->Uswap; /* update flow field */
    field->Uswap = exchanger; /* regain the used space as new space */

    LLy(field->Uswap, field->U, space, model, part, 0.5 * dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part,geometry);
    exchanger = field->U; /* preserve the address of U */
    field->U = field->Uswap; /* update flow field */
    field->Uswap = exchanger; /* regain the used space as new space */

    LLz(field->Uswap, field->U, space, model, part, 0.5 * dt);
    BoundaryCondtionsAndTreatments(field->Uswap, space, model, part, geometry);
    exchanger = field->U; /* preserve the address of U */
    field->U = field->Uswap; /* update flow field */
    field->Uswap = exchanger; /* regain the used space as new space */
    return 0;
}
static int LLz(Real *U, const Real *Un, const Space *space, const Model *model,
        const Partition *part, const Real dt)
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
                ComputeReconstructedFluxZ(Fhat, k, j, i, Un, space, model, dt);
                ComputeReconstructedFluxZ(Fhath, k - 1, j, i, Un, space, model, dt);
                ComputeViscousFluxGradientZ(gradG, k, j, i, Un, space, model);
                idx = idx * DIMU; /* change idx to field variable */
                for (int dim = 0; dim < DIMU; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int LLy(Real *U, const Real *Un, const Space *space, const Model *model,
        const Partition *part, const Real dt)
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
                ComputeReconstructedFluxY(Fhat, k, j, i, Un, space, model, dt);
                ComputeReconstructedFluxY(Fhath, k, j - 1, i, Un, space, model, dt);
                ComputeViscousFluxGradientY(gradG, k, j, i, Un, space, model);
                idx = idx * DIMU; /* change idx to field variable */
                for (int dim = 0; dim < DIMU; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int LLx(Real *U, const Real *Un, const Space *space, const Model *model,
        const Partition *part, const Real dt)
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
                ComputeReconstructedFluxX(Fhat, k, j, i, Un, space, model, dt);
                ComputeReconstructedFluxX(Fhath, k, j, i - 1, Un, space, model, dt);
                ComputeViscousFluxGradientX(gradG, k, j, i, Un, space, model);
                idx = idx * DIMU; /* change idx to field variable */
                for (int dim = 0; dim < DIMU; ++dim) {
                    U[idx+dim] = Un[idx+dim] - r * (Fhat[dim] - Fhath[dim]) + dt * gradG[dim];
                }
            }
        }
    }
    return 0;
}
static int ComputeReconstructedFluxZ(Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    TVDFluxZ(Fhat, k, j, i, U, space, model, dt);
    return 0;
}
static int ComputeReconstructedFluxY(Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    TVDFluxY(Fhat, k, j, i, U, space, model, dt);
    return 0;
}
static int ComputeReconstructedFluxX(Real Fhat[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model, const Real dt)
{
    TVDFluxX(Fhat, k, j, i, U, space, model, dt);
    return 0;
}
/*
 * Generally the viscous terms will only be discretized by central
 * difference scheme. However, for outermost layer of interior nodes, the
 * central scheme of them require the viscous flux vector at boundaries, 
 * because the corners of current computational domain haven't set with any
 * values, the viscous flux vector at boundaries can not be interpolated
 * with central scheme, therefore, need to change to forward or backward. 
 * This situation is the same to interior ghost nodes the central
 * difference scheme can't be applied because of lacking stencil. Thus,
 * they also need to be identified and modified.
 */
static int ComputeViscousFluxGradientZ(
        Real gradG[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real hG[DIMU] = {0.0}; /* viscous flux vector */
    Real Gh[DIMU] = {0.0}; /* viscous flux vector */
    Real h = 0; /* reciprocal of differencing distance */
    /* default is central scheme */
    int hl = k - 1;
    int hr = k + 1;
    if (space->ng == hl) { /* if left at boundary */
        ++hl;
    }
    if ((space->nz + space->ng - 1) == hr) { /* if right at boundary */
        --hr;
    }
    /* check ghost */
    const int idxl = IndexMath(hl, j, i, space);
    const int idxr = IndexMath(hr, j, i, space);
    if (OFFSET <= space->nodeFlag[idxl]) {
        ++hl;
    }
    if (OFFSET <= space->nodeFlag[idxr]) {
        --hr;
    }
    if (0 != (hr - hl)) { /* only do calculation when needed */
        h = space->ddz / (hr - hl);
        ComputeViscousFluxZ(hG, hl, j, i, U, space, model);
        ComputeViscousFluxZ(Gh, hr, j, i, U, space, model);
    }
    for (int row = 0; row < DIMU; ++row) {
        gradG[row] = h * (Gh[row] - hG[row]);
    }
    return 0;
}
static int ComputeViscousFluxGradientY(
        Real gradG[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real hG[DIMU] = {0.0}; /* viscous flux vector */
    Real Gh[DIMU] = {0.0}; /* viscous flux vector */
    Real h = 0; /* reciprocal of differencing distance */
    /* default is central scheme */
    int hl = j - 1;
    int hr = j + 1;
    if (space->ng == hl) { /* if left at boundary */
        ++hl;
    }
    if ((space->ny + space->ng - 1) == hr) { /* if right at boundary */
        --hr;
    }
    /* check ghost */
    const int idxl = IndexMath(k, hl, i, space);
    const int idxr = IndexMath(k, hr, i, space);
    if (OFFSET <= space->nodeFlag[idxl]) {
        ++hl;
    }
    if (OFFSET <= space->nodeFlag[idxr]) {
        --hr;
    }
    if (0 != (hr - hl)) { /* only do calculation when needed */
        h = space->ddy / (hr - hl);
        ComputeViscousFluxY(hG, k, hl, i, U, space, model);
        ComputeViscousFluxY(Gh, k, hr, i, U, space, model);
    }
    for (int row = 0; row < DIMU; ++row) {
        gradG[row] = h * (Gh[row] - hG[row]);
    }
    return 0;
}
static int ComputeViscousFluxGradientX(
        Real gradG[], const int k, const int j, const int i, 
        const Real *U, const Space *space, const Model *model)
{
    Real hG[DIMU] = {0.0}; /* viscous flux vector */
    Real Gh[DIMU] = {0.0}; /* viscous flux vector */
    Real h = 0; /* reciprocal of differencing distance */
    /* default is central scheme */
    int hl = i - 1;
    int hr = i + 1;
    if (space->ng == hl) { /* if left at boundary */
        ++hl;
    }
    if ((space->nx + space->ng - 1) == hr) { /* if right at boundary */
        --hr;
    }
    /* check ghost */
    const int idxl = IndexMath(k, j, hl, space);
    const int idxr = IndexMath(k, j, hr, space);
    if (OFFSET <= space->nodeFlag[idxl]) {
        ++hl;
    }
    if (OFFSET <= space->nodeFlag[idxr]) {
        --hr;
    }
    if (0 != (hr - hl)) { /* only do calculation when needed */
        h = space->ddx / (hr - hl);
        ComputeViscousFluxX(hG, k, j, hl, U, space, model);
        ComputeViscousFluxX(Gh, k, j, hr, U, space, model);
    }
    for (int row = 0; row < DIMU; ++row) {
        gradG[row] = h * (Gh[row] - hG[row]);
    }
    return 0;
}
static int ComputeViscousFluxZ(
        Real G[], const int k, const int j, const int i, 
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
static int ComputeViscousFluxY(
        Real G[], const int k, const int j, const int i, 
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
static int ComputeViscousFluxX(
        Real G[], const int k, const int j, const int i, 
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
/* a good practice: end file with a newline */

