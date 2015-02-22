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
    Real *rho = field->Un + 0 * space->kMax * space->jMax * space->iMax;
    Real *rho_u = field->Un + 1 * space->kMax * space->jMax * space->iMax;
    Real *rho_v = field->Un + 2 * space->kMax * space->jMax * space->iMax;
    Real *rho_w = field->Un + 3 * space->kMax * space->jMax * space->iMax;
    Real *rho_eT = field->Un + 4 * space->kMax * space->jMax * space->iMax;
    /*
     * Decompose the primitive field variable into each component.
     */
    Real *rhoPri = field->Uo + 0 * space->kMax * space->jMax * space->iMax;
    Real *u = field->Uo + 1 * space->kMax * space->jMax * space->iMax;
    Real *v = field->Uo + 2 * space->kMax * space->jMax * space->iMax;
    Real *w = field->Uo + 3 * space->kMax * space->jMax * space->iMax;
    Real *p = field->Uo + 4 * space->kMax * space->jMax * space->iMax;
    Real *T = field->Uo + 5 * space->kMax * space->jMax * space->iMax;
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
                rhoPri[idx] = rho[idx];
                u[idx] = rho_u[idx] / rho[idx];
                v[idx] = rho_v[idx] / rho[idx];
                w[idx] = rho_w[idx] / rho[idx];
                p[idx] = (flow->gamma - 1) * (rho_eT[idx] - 0.5 * rho[idx] * 
                        (u[idx] * u[idx] + v[idx] * v[idx] + w[idx] * w[idx]));
                T[idx] = p[idx] / (rho[idx] * flow->gasR);
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
    Real *rho = field->Uo + 0 * space->kMax * space->jMax * space->iMax;
    Real *u = field->Uo + 1 * space->kMax * space->jMax * space->iMax;
    Real *v = field->Uo + 2 * space->kMax * space->jMax * space->iMax;
    Real *w = field->Uo + 3 * space->kMax * space->jMax * space->iMax;
    Real *p = field->Uo + 4 * space->kMax * space->jMax * space->iMax;
    /*
     * Auxiliary variables
     */
    Real velocity = 0;
    Real velocityMax = 1e-100;
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
        return 1e100;
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

