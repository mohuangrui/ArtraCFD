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
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_LINEAR_SYSTEM_H_ /* if this is the first definition */
#define ARTRACFD_LINEAR_SYSTEM_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Linear System Solver in Matrix Form
 *
 * Function
 *      General matrix linear system solver for AX = B. 
 *      Dimension of A is n x n; dimension of X and B is n x m;
 *      A is the system matrix, X is the solution matrix with each column vector
 *      stands for a solution vector of a single linear system A * x = b, B is 
 *      the right hand matrix with each column vector stands for a right hand 
 *      vector of a single linear system. If B is a unit matrix, then solution
 *      matrix X is the inverse matrix of A. In practice, right hand matrix 
 *      and solution matrix can share the same storage even they are declared
 *      as restricted pointers here. The matrix A will be transformed into
 *      its LU decomposition.
 */
extern int MatrixLinearSystemSolver(const int n, Real A[restrict][n],
        const int m, Real X[restrict][m], Real B[restrict][m]);
#endif
/* a good practice: end file with a newline */

