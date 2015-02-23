/****************************************************************************
 * Some common functions for CFD                                            *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines some common functions for CFD.                       *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "cfdcommons.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ComputePrimitiveByConservative(Field *field, const Space *space, const Flow *flow)
{
    /*
     * Decompose the conservative field variable into each component.
     */
    const Real *Un[5] = {
        field->Un + 0 * space->nMax,
        field->Un + 1 * space->nMax,
        field->Un + 2 * space->nMax,
        field->Un + 3 * space->nMax,
        field->Un + 4 * space->nMax};
    /*
     * Decompose the primitive field variable into each component.
     */
    Real *rho = field->Uo + 0 * space->nMax;
    Real *u = field->Uo + 1 * space->nMax;
    Real *v = field->Uo + 2 * space->nMax;
    Real *w = field->Uo + 3 * space->nMax;
    Real *p = field->Uo + 4 * space->nMax;
    Real *T = field->Uo + 5 * space->nMax;
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* calculated index */
    for (k = 0; k < space->kMax; ++k) {
        for (j = 0; j < space->jMax; ++j) {
            for (i = 0; i < space->iMax; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (space->ghostFlag[idx] == -1) { /* if it's solid node */
                    continue;
                }
                rho[idx] = Un[0][idx];
                u[idx] = Un[1][idx] / rho[idx];
                v[idx] = Un[2][idx] / rho[idx];
                w[idx] = Un[3][idx] / rho[idx];
                p[idx] = (flow->gamma - 1) * (Un[4][idx] - 0.5 * rho[idx] * 
                        (u[idx] * u[idx] + v[idx] * v[idx] + w[idx] * w[idx]));
                T[idx] = p[idx] / (rho[idx] * flow->gasR);
            }
        }
    }
    return 0;
}
int ComputeConservativeByPrimitive(Field *field, const Space *space, const Flow *flow)
{
    /*
     * Decompose the conservative field variable into each component.
     */
    Real *Un[5] = {
        field->Un + 0 * space->nMax,
        field->Un + 1 * space->nMax,
        field->Un + 2 * space->nMax,
        field->Un + 3 * space->nMax,
        field->Un + 4 * space->nMax};
    /*
     * Decompose the primitive field variable into each component.
     */
    const Real *rho = field->Uo + 0 * space->nMax;
    const Real *u = field->Uo + 1 * space->nMax;
    const Real *v = field->Uo + 2 * space->nMax;
    const Real *w = field->Uo + 3 * space->nMax;
    const Real *p = field->Uo + 4 * space->nMax;
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* calculated index */
    for (k = 0; k < space->kMax; ++k) {
        for (j = 0; j < space->jMax; ++j) {
            for (i = 0; i < space->iMax; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (space->ghostFlag[idx] == -1) { /* if it's solid node */
                    continue;
                }
                Un[0][idx] = rho[idx];
                Un[1][idx] = rho[idx] * u[idx];
                Un[2][idx] = rho[idx] * v[idx];
                Un[3][idx] = rho[idx] * w[idx];
                Un[4][idx] = p[idx] / (flow->gamma - 1) + 
                    0.5 * rho[idx] * (u[idx] * u[idx] + v[idx] * v[idx] + w[idx] * w[idx]);
            }
        }
    }
    return 0;
}
int ComputeNonviscousFlux(const Field *field, Flux *flux, const Space *space)
{
    /*
     * Decompose the nonviscous flux variables into each component
     */
    Real *Fx[5] = {
        flux->Fx + 0 * space->nMax,
        flux->Fx + 1 * space->nMax,
        flux->Fx + 2 * space->nMax,
        flux->Fx + 3 * space->nMax,
        flux->Fx + 4 * space->nMax};
    Real *Fy[5] = {
        flux->Fy + 0 * space->nMax,
        flux->Fy + 1 * space->nMax,
        flux->Fy + 2 * space->nMax,
        flux->Fy + 3 * space->nMax,
        flux->Fy + 4 * space->nMax};
    Real *Fz[5] = {
        flux->Fz + 0 * space->nMax,
        flux->Fz + 1 * space->nMax,
        flux->Fz + 2 * space->nMax,
        flux->Fz + 3 * space->nMax,
        flux->Fz + 4 * space->nMax};
    /*
     * Decompose the primitive field variable into each component.
     */
    const Real *rho = field->Uo + 0 * space->nMax;
    const Real *u = field->Uo + 1 * space->nMax;
    const Real *v = field->Uo + 2 * space->nMax;
    const Real *w = field->Uo + 3 * space->nMax;
    const Real *p = field->Uo + 4 * space->nMax;
    const Real *rho_eT = field->Un + 4 * space->nMax;
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* calculated index */
    for (k = 0; k < space->kMax; ++k) {
        for (j = 0; j < space->jMax; ++j) {
            for (i = 0; i < space->iMax; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (space->ghostFlag[idx] == -1) { /* if it's solid node */
                    continue;
                }
                Fx[0][idx] = rho[idx] * u[idx];
                Fx[1][idx] = rho[idx] * u[idx] * u[idx] + p[idx];
                Fx[2][idx] = rho[idx] * u[idx] * v[idx];
                Fx[3][idx] = rho[idx] * u[idx] * w[idx];
                Fx[4][idx] = (rho_eT[idx] + p[idx]) * u[idx];

                Fy[0][idx] = rho[idx] * v[idx];
                Fy[1][idx] = rho[idx] * v[idx] * u[idx];
                Fy[2][idx] = rho[idx] * v[idx] * v[idx] + p[idx];
                Fy[3][idx] = rho[idx] * v[idx] * w[idx];
                Fy[4][idx] = (rho_eT[idx] + p[idx]) * v[idx];

                Fz[0][idx] = rho[idx] * w[idx];
                Fz[1][idx] = rho[idx] * w[idx] * u[idx];
                Fz[2][idx] = rho[idx] * w[idx] * v[idx];
                Fz[3][idx] = rho[idx] * w[idx] * w[idx] + p[idx];
                Fz[4][idx] = (rho_eT[idx] + p[idx]) * w[idx];
            }
        }
    }
    return 0;
}
int ComputeViscousFlux(const Field *field, Flux *flux, const Space *space, const Flow *flow)
{
    /*
     * Decompose the viscous flux variables into each component
     */
    Real *Gx[5] = {
        flux->Gx + 0 * space->nMax,
        flux->Gx + 1 * space->nMax,
        flux->Gx + 2 * space->nMax,
        flux->Gx + 3 * space->nMax,
        flux->Gx + 4 * space->nMax};
    Real *Gy[5] = {
        flux->Gy + 0 * space->nMax,
        flux->Gy + 1 * space->nMax,
        flux->Gy + 2 * space->nMax,
        flux->Gy + 3 * space->nMax,
        flux->Gy + 4 * space->nMax};
    Real *Gz[5] = {
        flux->Gz + 0 * space->nMax,
        flux->Gz + 1 * space->nMax,
        flux->Gz + 2 * space->nMax,
        flux->Gz + 3 * space->nMax,
        flux->Gz + 4 * space->nMax};
    /*
     * Decompose the primitive field variable into each component.
     */
    const Real *u = field->Uo + 1 * space->nMax;
    const Real *v = field->Uo + 2 * space->nMax;
    const Real *w = field->Uo + 3 * space->nMax;
    const Real *T = field->Uo + 5 * space->nMax;
    /*
     * Auxiliary variables
     */
    Real divV = 0; /* divergence of velocity V */
    Real du_dx = 0; /* partial u partial x */
    Real du_dy = 0;
    Real du_dz = 0;
    Real dv_dx = 0; /* partial v partial x */
    Real dv_dy = 0;
    Real dv_dz = 0;
    Real dw_dx = 0; /* partial w partial x */
    Real dw_dy = 0;
    Real dw_dz = 0;
    Real dT_dx = 0; /* partial T partial x */
    Real dT_dy = 0;
    Real dT_dz = 0;
    const Real dx = MinPositive(space->dx, -1);
    const Real dy = MinPositive(space->dy, -1);
    const Real dz = MinPositive(space->dz, -1);
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* calculated index */
    int idxW = 0; /* index at West */
    int idxE = 0; /* index at East */
    int idxS = 0; /* index at South */
    int idxN = 0; /* index at North */
    int idxF = 0; /* index at Front */
    int idxB = 0; /* index at Back */
    for (k = 0; k < space->kMax; ++k) {
        for (j = 0; j < space->jMax; ++j) {
            for (i = 0; i < space->iMax; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (space->ghostFlag[idx] == -1) { /* if it's solid node */
                    continue;
                }

                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxE = (k * space->jMax + j) * space->iMax + i + 1;
                idxS = (k * space->jMax + j - 1) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                idxF = ((k - 1) * space->jMax + j) * space->iMax + i;
                idxB = ((k + 1) * space->jMax + j) * space->iMax + i;

                if (k == 0) {
                    du_dz = (u[idxB] - u[idx]) / dz;
                    dv_dz = (v[idxB] - v[idx]) / dz;
                    dw_dz = (w[idxB] - w[idx]) / dz;
                    dT_dz = (T[idxB] - T[idx]) / dz;
                } else {
                    if (k == space->kMax - 1) {
                        du_dz = (u[idx] - u[idxF]) / dz;
                        dv_dz = (v[idx] - v[idxF]) / dz;
                        dw_dz = (w[idx] - w[idxF]) / dz;
                        dT_dz = (T[idx] - T[idxF]) / dz;
                    } else {
                        du_dz = (u[idxB] - u[idxF]) / (2 * dz);
                        dv_dz = (v[idxB] - v[idxF]) / (2 * dz);
                        dw_dz = (w[idxB] - w[idxF]) / (2 * dz);
                        dT_dz = (T[idxB] - T[idxF]) / (2 * dz);
                    }
                }
                if (j == 0) {
                    du_dy = (u[idxN] - u[idx]) / dy;
                    dv_dy = (v[idxN] - v[idx]) / dy;
                    dw_dy = (w[idxN] - w[idx]) / dy;
                    dT_dy = (T[idxN] - T[idx]) / dy;
                } else {
                    if (j == space->jMax - 1) {
                        du_dy = (u[idx] - u[idxS]) / dy;
                        dv_dy = (v[idx] - v[idxS]) / dy;
                        dw_dy = (w[idx] - w[idxS]) / dy;
                        dT_dy = (T[idx] - T[idxS]) / dy;
                    } else {
                        du_dy = (u[idxN] - u[idxS]) / (2 * dy);
                        dv_dy = (v[idxN] - v[idxS]) / (2 * dy);
                        dw_dy = (w[idxN] - w[idxS]) / (2 * dy);
                        dT_dy = (T[idxN] - T[idxS]) / (2 * dy);
                    }
                }
                if (i == 0) {
                    du_dx = (u[idxE] - u[idx]) / dx;
                    dv_dx = (v[idxE] - v[idx]) / dx;
                    dw_dx = (w[idxE] - w[idx]) / dx;
                    dT_dx = (T[idxE] - T[idx]) / dx;
                } else {
                    if (i == space->iMax - 1) {
                        du_dx = (u[idx] - u[idxW]) / dx;
                        dv_dx = (v[idx] - v[idxW]) / dx;
                        dw_dx = (w[idx] - w[idxW]) / dx;
                        dT_dx = (T[idx] - T[idxW]) / dx;
                    } else {
                        du_dx = (u[idxE] - u[idxW]) / (2 * dx);
                        dv_dx = (v[idxE] - v[idxW]) / (2 * dx);
                        dw_dx = (w[idxE] - w[idxW]) / (2 * dx);
                        dT_dx = (T[idxE] - T[idxW]) / (2 * dx);
                    }
                }

                divV = du_dx + dv_dy + dw_dz;

                Gx[0][idx] = 0;
                Gx[1][idx] = flow->mu * (2 * du_dx - (2/3) * divV);
                Gx[2][idx] = flow->mu * (du_dy + dv_dx);
                Gx[3][idx] = flow->mu * (du_dz + dw_dx);
                Gx[4][idx] = flow->heatK * dT_dx + 
                    u[idx] * Gx[1][idx] + v[idx] * Gx[2][idx] + w[idx] * Gx[3][idx];

                Gy[0][idx] = 0;
                Gy[1][idx] = Gx[2][idx];
                Gy[2][idx] = flow->mu * (2 * dv_dy - (2/3) * divV);
                Gy[3][idx] = flow->mu * (dv_dz + dw_dy);
                Gy[4][idx] = flow->heatK * dT_dy + 
                    u[idx] * Gy[1][idx] + v[idx] * Gy[2][idx] + w[idx] * Gy[3][idx];

                Gz[0][idx] = 0;
                Gz[1][idx] = Gx[3][idx];
                Gz[2][idx] = Gy[3][idx];
                Gz[3][idx] = flow->mu * (2 * dw_dz - (2/3) * divV);
                Gz[4][idx] = flow->heatK * dT_dz + 
                    u[idx] * Gz[1][idx] + v[idx] * Gz[2][idx] + w[idx] * Gz[3][idx];
            }
        }
    }
    return 0;
}
Real ComputeTimeStepByCFL(Field *field, const Space *space, const Time *time, 
        const Partition *part, const Flow *flow)
{
    /*
     * Decompose the primitive field variable into each component.
     */
    Real *rho = field->Uo + 0 * space->nMax;
    Real *u = field->Uo + 1 * space->nMax;
    Real *v = field->Uo + 2 * space->nMax;
    Real *w = field->Uo + 3 * space->nMax;
    Real *p = field->Uo + 4 * space->nMax;
    /*
     * Auxiliary variables
     */
    Real velocity = 0;
    Real velocityMax = 1e-38;
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* calculated index */
    for (k = part->kSub[12]; k < part->kSup[12]; ++k) {
        for (j = part->jSub[12]; j < part->jSup[12]; ++j) {
            for (i = part->iSub[12]; i < part->iSup[12]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (space->ghostFlag[idx] == -1) { /* if it's solid node */
                    continue;
                }
                velocity = sqrt(flow->gamma * p[idx] / rho[idx]) + 
                    sqrt(u[idx] * u[idx] + v[idx] * v[idx] + w[idx] *w[idx]);
                if (velocityMax < velocity) {
                    velocityMax = velocity;
                }
            }
        }
    }
    return time->numCFL * MinPositive(space->dx, MinPositive(space->dy, space->dz)) / velocityMax;
}
Real MinPositive(Real valueA, Real valueB)
{
    if ((valueA <= 0) && (valueB <= 0)) {
        return 1e38;
    }
    if (valueA <= 0) {
        return valueB;
    }
    if (valueB <= 0) {
        return valueA;
    }
    return Min(valueA, valueB);
}
Real Min(Real valueA, Real valueB)
{
    if (valueA < valueB) {
        return valueA;
    }
    return valueB;
}
/* a good practice: end file with a newline */

