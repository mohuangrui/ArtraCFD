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
#include "weno.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    R = 3, /* WENO r */
    CN = 2, /* position index of the center node in stencil */
} WENOConst;
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static Real Square(const Real);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Jiang, G.S. and Shu, C.W., 1996. Efficient Implementation of Weighted
 * ENO Schemes. Journal of Computational Physics, 126(1), pp.202-228.
 */
void WENO5(Real F[restrict][DIMU], Real Fhat[restrict])
{
    Real omega[R]; /* weights */
    Real q[R]; /* q vectors */
    Real IS[R]; /* smoothness measurements */
    Real alpha[R];
    const Real C[R] = {1.0 / 10.0, 6.0 / 10.0, 3.0 / 10.0};
    const Real epsilon = 1.0e-6;
    for (int r = 0; r < DIMU; ++r) {
        IS[0] = (13.0 / 12.0) * Square(F[CN-2][r] - 2.0 * F[CN-1][r] + F[CN][r]) +
            (1.0 / 4.0) * Square(F[CN-2][r] - 4.0 * F[CN-1][r] + 3.0 * F[CN][r]);
        IS[1] = (13.0 / 12.0) * Square(F[CN-1][r] - 2.0 * F[CN][r] + F[CN+1][r]) +
            (1.0 / 4.0) * Square(F[CN-1][r] - F[CN+1][r]);
        IS[2] = (13.0 / 12.0) * Square(F[CN][r] - 2.0 * F[CN+1][r] + F[CN+2][r]) +
            (1.0 / 4.0) * Square(3.0 * F[CN][r] - 4.0 * F[CN+1][r] + F[CN+2][r]);
        alpha[0] = C[0] / Square(epsilon + IS[0]);
        alpha[1] = C[1] / Square(epsilon + IS[1]);
        alpha[2] = C[2] / Square(epsilon + IS[2]);
        omega[0] = alpha[0] / (alpha[0] + alpha[1] + alpha[2]);
        omega[1] = alpha[1] / (alpha[0] + alpha[1] + alpha[2]);
        omega[2] = alpha[2] / (alpha[0] + alpha[1] + alpha[2]);
        q[0] = (1.0 / 6.0) * (2.0 * F[CN-2][r] - 7.0 * F[CN-1][r] + 11.0 * F[CN][r]);
        q[1] = (1.0 / 6.0) * (-F[CN-1][r] + 5.0 * F[CN][r] + 2.0 * F[CN+1][r]);
        q[2] = (1.0 / 6.0) * (2.0 * F[CN][r] + 5.0 * F[CN+1][r] - F[CN+2][r]);
        Fhat[r] = omega[0] * q[0] + omega[1] * q[1] + omega[2] * q[2];
    }
    return;
}
static Real Square(const Real x)
{
    return x * x;
}
/* a good practice: end file with a newline */

