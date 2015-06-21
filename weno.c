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
 * Function Pointers
 ****************************************************************************/
/*
 * Function pointers are useful for implementing a form of polymorphism.
 * They are mainly used to reduce or avoid switch statement. Pointers to
 * functions can get rather messy. Declaring a typedel to a function pointer
 * generally clarifies the code.
 */
typedef int (*ConvectiveFluxComputer)(Real [], const int, const int, const int, 
        const Real *, const Space *, const Model *);
typedef int (*EigenvectorSpaceRComputer)(Real [][DIMU], const int, const int,
        const int, const Real *, const Space *, const Model *);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WENO(const int s, Real Fhat[], const Real r, const int k, const int j,
        const int i, const Real *U, const Space *space, const Model *model)
{
    Real Fplus[DIMU] = {0.0}; /* forward numerical flux */
    Real Fminus[DIMU] = {0.0}; /* backward numerical flux */
    for (int row = 0; row < DIMU; ++row) {
        Fhat[row] = Fplus[row] + Fminus[row];
    }
    return 0;
}
/* a good practice: end file with a newline */

