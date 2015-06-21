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
#include "weno.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    R = 3, /* WENO r */
    N = 2, /* N = (r - 1) */
    CEN = 2, /* center index of stencil */
    NSTENCIL = 5, /* number of nodes in a stencil = (2r - 1) */
} WENOConstants;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NumericalFlux(Real [], Real [][DIMU]);
static int Omega(Real [][DIMU], Real [][DIMU]);
static int Q(Real [][DIMU], Real [][DIMU]);
static int ComputeIS(Real [][DIMU], Real [][DIMU]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WENO(const int s, Real Fhat[], const Real r, const int k, const int j,
        const int i, const Real *U, const Space *space, const Model *model)
{
    Real FhatPlus[DIMU] = {0.0}; /* forward numerical flux */
    Real FhatMinus[DIMU] = {0.0}; /* backward numerical flux */
    Real Fplus[NSTENCIL][DIMU] = {{0.0}}; /* forward flux of each stencil */
    Real Fminus[NSTENCIL][DIMU] = {{0.0}}; /* backward flux of each stencil */
    const int h[DIMS][DIMS] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* direction indicator */
    for (int n = -N; n <= N; ++n) {
        /* forward and backward fluxes should be symmetrically organized */
        FluxVectorSplitting(s, Fplus[CEN+n], Fminus[CEN-n], k + h[s][Z] * n,
                j + h[s][Y] * n, i + h[s][X] * n, U, space, model);
    }
    NumericalFlux(FhatPlus, Fplus);
    NumericalFlux(FhatMinus, Fminus);
    for (int row = 0; row < DIMU; ++row) {
        Fhat[row] = FhatPlus[row] + FhatMinus[row];
    }
    return (!r); /* to avoid unsed parameter warning */
}
static int NumericalFlux(Real Fhat[], Real F[][DIMU])
{
    Real omega[R][DIMU] = {{0.0}}; /* weights */
    Real q[R][DIMU] = {{0.0}}; /* q vectors */
    Omega(omega, F);
    Q(q, F);
    for (int row = 0; row < DIMU; ++row) {
        Fhat[row] = omega[0][row] * q[0][row] + omega[1][row] * q[1][row] + omega[2][row] * q[2][row];
    }
    return 0;
}
static int Omega(Real omega[][DIMU], Real F[][DIMU])
{
    const Real C[R] = {0.1, 0.6, 0.3};
    const Real epsilon = 1.0e-6;
    Real IS[R][DIMU] = {{0.0}};
    Real alpha[R][DIMU] = {{0.0}};
    ComputeIS(IS, F);
    for (int row = 0; row < DIMU; ++row) {
        alpha[0][row] = C[0] / pow((epsilon + IS[0][row]), 2.0);
        alpha[1][row] = C[1] / pow((epsilon + IS[1][row]), 2.0);
        alpha[2][row] = C[2] / pow((epsilon + IS[2][row]), 2.0);
    }
    for (int row = 0; row < DIMU; ++row) {
        omega[0][row] = alpha[0][row] / (alpha[0][row] + alpha[1][row] + alpha[2][row]);
        omega[1][row] = alpha[1][row] / (alpha[0][row] + alpha[1][row] + alpha[2][row]);
        omega[2][row] = alpha[2][row] / (alpha[0][row] + alpha[1][row] + alpha[2][row]);
    }
    return 0;
}
static int ComputeIS(Real IS[][DIMU], Real F[][DIMU])
{
    for (int row = 0; row < DIMU; ++row) {
        IS[0][row] = (13.0 / 12.0) * pow((F[CEN-2][row] - 2.0 * F[CEN-1][row] + F[CEN][row]), 2.0) + 
            (1.0 / 4.0) * pow((F[CEN-2][row] - 4.0 * F[CEN-1][row] + 3.0 * F[CEN][row]), 2.0);
        IS[1][row] = (13.0 / 12.0) * pow((F[CEN-1][row] - 2.0 * F[CEN][row] + F[CEN+1][row]), 2.0) +
            (1.0 / 4.0) * pow((F[CEN-1][row] - F[CEN+1][row]), 2.0);
        IS[2][row] = (13.0 / 12.0) * pow((F[CEN][row] - 2.0 * F[CEN+1][row] + F[CEN+2][row]), 2.0) +
            (1.0 / 4.0) * pow((3.0 * F[CEN][row] - 4.0 * F[CEN+1][row] + F[CEN+2][row]), 2.0);
    }
    return 0;
}
static int Q(Real q[][DIMU], Real F[][DIMU])
{
    for (int row = 0; row < DIMU; ++row) {
        q[0][row] = (2.0 * F[CEN-2][row] - 7.0 * F[CEN-1][row] + 11.0 * F[CEN][row]) / 6.0;
        q[1][row] = (-F[CEN-1][row] + 5.0 * F[CEN][row] + 2.0 * F[CEN+1][row]) / 6.0;
        q[2][row] = (2.0 * F[CEN][row] + 5.0 * F[CEN+1][row] - F[CEN+2][row]) / 6.0;
    }
    return 0;
}
/* a good practice: end file with a newline */

