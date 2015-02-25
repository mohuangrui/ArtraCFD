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
int ComputeNonviscousFlux(const Field *field, Flux *flux, const Space *space, const Flow *flow)
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
     * Define the primitive field variables.
     */
    Real rho = 0; 
    Real u = 0;
    Real v = 0;
    Real w = 0;
    Real p = 0;
    Real eT = 0;
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
                rho = Un[0][idx];
                u = Un[1][idx] / rho;
                v = Un[2][idx] / rho;
                w = Un[3][idx] / rho;
                eT = Un[4][idx] / rho;
                p = (flow->gamma - 1) * rho * (eT - 0.5 * (u * u + v * v + w * w));

                Fx[0][idx] = rho * u;
                Fx[1][idx] = rho * u * u + p;
                Fx[2][idx] = rho * u * v;
                Fx[3][idx] = rho * u * w;
                Fx[4][idx] = (rho * eT + p) * u;

                Fy[0][idx] = rho * v;
                Fy[1][idx] = rho * v * u;
                Fy[2][idx] = rho * v * v + p;
                Fy[3][idx] = rho * v * w;
                Fy[4][idx] = (rho * eT + p) * v;

                Fz[0][idx] = rho * w;
                Fz[1][idx] = rho * w * u;
                Fz[2][idx] = rho * w * v;
                Fz[3][idx] = rho * w * w + p;
                Fz[4][idx] = (rho * eT + p) * w;
            }
        }
    }
    return 0;
}
int ComputeViscousFlux(const Field *field, Flux *flux, const Space *space, const Flow *flow)
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
     * Define the primitive field variables.
     */
    Real rho = 0; 
    Real rho_h = 0; 
    Real u = 0;
    Real u_h = 0;
    Real v = 0;
    Real v_h = 0;
    Real w = 0;
    Real w_h = 0;
    Real eT = 0;
    Real eT_h = 0;
    Real T = 0;
    Real T_h = 0;
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
    const Real dx = MinPositive(space->dx, -1); /* needed when use as denominator */
    const Real dy = MinPositive(space->dy, -1); /* needed when use as denominator */
    const Real dz = MinPositive(space->dz, -1); /* needed when use as denominator */
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
    /*
     * Generally the viscous terms will only be discretized by central
     * difference scheme, therefore, only the viscous variables at 
     * the boudary region and interior region need to be calculated, 
     * there is no need to do this for outer ghost cells.
     */
    for (k = space->ng; k < space->kMax - space->ng; ++k) {
        for (j = space->ng; j < space->jMax - space->ng; ++j) {
            for (i = space->ng; i < space->iMax - space->ng; ++i) {
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

                /* calculate derivatives in z direction */
                rho_h = Un[0][idxB];
                u_h = Un[1][idxB] / rho_h;
                v_h = Un[2][idxB] / rho_h;
                w_h = Un[3][idxB] / rho_h;
                eT_h = Un[4][idxB] / rho_h;
                T_h = (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h)) / flow->cv;

                rho = Un[0][idxF];
                u = Un[1][idxF] / rho;
                v = Un[2][idxF] / rho;
                w = Un[3][idxF] / rho;
                eT = Un[4][idxF] / rho;
                T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

                du_dz = (u_h - u) / (2 * dz);
                dv_dz = (v_h - v) / (2 * dz);
                dw_dz = (w_h - w) / (2 * dz);
                dT_dz = (T_h - T) / (2 * dz);

                /* calculate derivatives in y direction */
                rho_h = Un[0][idxN];
                u_h = Un[1][idxN] / rho_h;
                v_h = Un[2][idxN] / rho_h;
                w_h = Un[3][idxN] / rho_h;
                eT_h = Un[4][idxN] / rho_h;
                T_h = (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h)) / flow->cv;

                rho = Un[0][idxS];
                u = Un[1][idxS] / rho;
                v = Un[2][idxS] / rho;
                w = Un[3][idxS] / rho;
                eT = Un[4][idxS] / rho;
                T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

                du_dy = (u_h - u) / (2 * dy);
                dv_dy = (v_h - v) / (2 * dy);
                dw_dy = (w_h - w) / (2 * dy);
                dT_dy = (T_h - T) / (2 * dy);

                /* calculate derivatives in x direction */
                rho_h = Un[0][idxE];
                u_h = Un[1][idxE] / rho_h;
                v_h = Un[2][idxE] / rho_h;
                w_h = Un[3][idxE] / rho_h;
                eT_h = Un[4][idxE] / rho_h;
                T_h = (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h)) / flow->cv;

                rho = Un[0][idxW];
                u = Un[1][idxW] / rho;
                v = Un[2][idxW] / rho;
                w = Un[3][idxW] / rho;
                eT = Un[4][idxW] / rho;
                T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;

                du_dx = (u_h - u) / (2 * dx);
                dv_dx = (v_h - v) / (2 * dx);
                dw_dx = (w_h - w) / (2 * dx);
                dT_dx = (T_h - T) / (2 * dx);

                /* regain the primitive variables in current point */
                rho = Un[0][idx];
                u = Un[1][idx] / rho;
                v = Un[2][idx] / rho;
                w = Un[3][idx] / rho;

                divV = du_dx + dv_dy + dw_dz;

                Gx[0][idx] = 0;
                Gx[1][idx] = flow->mu * (2 * du_dx - (2/3) * divV);
                Gx[2][idx] = flow->mu * (du_dy + dv_dx);
                Gx[3][idx] = flow->mu * (du_dz + dw_dx);
                Gx[4][idx] = flow->heatK * dT_dx + 
                    u * Gx[1][idx] + v * Gx[2][idx] + w * Gx[3][idx];

                Gy[0][idx] = 0;
                Gy[1][idx] = Gx[2][idx];
                Gy[2][idx] = flow->mu * (2 * dv_dy - (2/3) * divV);
                Gy[3][idx] = flow->mu * (dv_dz + dw_dy);
                Gy[4][idx] = flow->heatK * dT_dy + 
                    u * Gy[1][idx] + v * Gy[2][idx] + w * Gy[3][idx];

                Gz[0][idx] = 0;
                Gz[1][idx] = Gx[3][idx];
                Gz[2][idx] = Gy[3][idx];
                Gz[3][idx] = flow->mu * (2 * dw_dz - (2/3) * divV);
                Gz[4][idx] = flow->heatK * dT_dz + 
                    u * Gz[1][idx] + v * Gz[2][idx] + w * Gz[3][idx];
            }
        }
    }
    return 0;
}
Real ComputeTimeStepByCFL(const Field *field, const Space *space, const Time *time, 
        const Partition *part, const Flow *flow)
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
     * Define the primitive field variables.
     */
    Real rho = 0; 
    Real u = 0;
    Real v = 0;
    Real w = 0;
    Real p = 0;
    Real eT = 0;
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
                rho = Un[0][idx];
                u = Un[1][idx] / rho;
                v = Un[2][idx] / rho;
                w = Un[3][idx] / rho;
                eT = Un[4][idx] / rho;
                p = (flow->gamma - 1) * rho * (eT - 0.5 * (u * u + v * v + w * w));

                velocity = sqrt(flow->gamma * p / rho) + 
                    Max(u, Max(v, w));
                if (velocityMax < velocity) {
                    velocityMax = velocity;
                }
            }
        }
    }
    return time->numCFL * MinPositive(space->dx, MinPositive(space->dy, space->dz)) / velocityMax;
}
Real MinPositive(const Real valueA, const Real valueB)
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
Real Min(const Real valueA, const Real valueB)
{
    if (valueA < valueB) {
        return valueA;
    }
    return valueB;
}
Real Max(const Real valueA, const Real valueB)
{
    if (valueA > valueB) {
        return valueA;
    }
    return valueB;
}
/* a good practice: end file with a newline */

