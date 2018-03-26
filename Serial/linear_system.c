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
#include "linear_system.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int LUFactorization(const int n, Real A[restrict][n], int permute[restrict]);
static int FactorizedLinearSystemSolver(const int n, Real L[restrict][n], Real U[restrict][n],
        Real x[], Real b[], const int permute[restrict]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int MatrixLinearSystemSolver(const int n, Real A[restrict][n],
        const int m, Real X[restrict][m], Real B[restrict][m])
{
    int permute[n]; /* record the permutation information */
    Real rhs[n]; /* transfer data into column vector */
    LUFactorization(n, A, permute);
    for (int col = 0; col < m; ++col) { /* solve column by column */
        for (int row = 0; row < n; ++row) {
            rhs[row] = B[row][col]; /* obtain each right hand side vector */
        }
        FactorizedLinearSystemSolver(n, A, A, rhs, rhs, permute);
        for (int row = 0; row < n; ++row) { /* save solution vector */
            X[row][col] = rhs[row];
        }
    }
    return 0;
}
/*
 * Perform A = LU factorization, A is a n x n matrix.
 * After factorization, both L and U are stored into the storage space of A,
 * The returned A has its upper triangle as U, and the part beneath the 
 * diagonal is equal to L. The missing diagonal elements of L are all 1.
 * Crout's algorithm is employed for factorization.
 *
 * Partial pivoting is used to stabilize the algorithm, that is, only row-wise
 * permutation will be employed. The permutations are recorded in the integer 
 * vector "permute", which will be used to reorder the right hand side vectors
 * before solving the factorized linear system.
 *
 * Press, William H. Numerical recipes: The art of scientific computing. 
 * Cambridge university press.
 */
static int LUFactorization(const int n, Real A[restrict][n], int permute[restrict])
{
    const Real epsilon = 1.0e-15; /* a small number for singularity check */
    const Real zero = 0.0;
    const Real one = 1.0;
    Real temp = 0.0; /* auxiliary variable */
    Real maximum = 0.0; /* store the maximum value in column */
    Real scale[n]; /* stores the implicit scaling of each row with variable length array */
    int rowMax = 0; /* record the row number of current pivot element */
    int sign = 1; /* used for determine the sign of determinant after pivoting */
    /*
     * Loop over rows to get the implicit scaling information.
     */
    for (int row = 0; row < n; ++row) {
        maximum = zero;
        for (int col = 0; col < n; ++col) {
            temp = fabs(A[row][col]);
            if (temp > maximum) {
                maximum = temp;
            }
        }
        if (zero == maximum) {
            FatalError("singular matrix in LU factorization...");
        }
        scale[row] = one / maximum; /* save the scaling */
    }
    /*
     * Do LU factorization with partial pivoting.
     */
    for (int loop = 0; loop < n; ++loop) {
        maximum = zero; /* initialize for the search of largest pivot element */
        rowMax = loop; /* initialize the pivot position to current row */
        for (int row = loop; row < n; ++row) { /* search pivot element for current loop */
            temp = scale[row] * fabs(A[row][loop]);
            if (temp > maximum) {
                maximum = temp;
                rowMax = row;
            }
        }
        if (loop != rowMax) { /* interchange rows if required */
            for (int col = 0; col < n; ++col) {
                temp = A[rowMax][col];
                A[rowMax][col] = A[loop][col];
                A[loop][col] = temp;
            }
            sign = -sign; /* change the parity of sign */
            scale[rowMax] = scale[loop]; /* replace the scale factor */
        }
        permute[loop] = rowMax; /* record the permutation */
        if (zero == A[loop][loop]) { /* substitute zero pivot element with epsilon */
            A[loop][loop] = epsilon;
        }
        for (int row = loop + 1; row < n; ++row) {
            A[row][loop] = A[row][loop] / A[loop][loop]; /* divide column element with pivot element */
            for (int col = loop + 1; col < n; ++col) { /* reduce remaining submatrix */
                A[row][col] = A[row][col] - A[row][loop] * A[loop][col];
            }
        }
    }
    return sign;
}
/*
 * Solve a factorized linear system LU * x = b. x is the solution vector, b is
 * the right hand side vector. Vector x and b can refer to the same storage
 * space, in this situation, solution vector overwrites the right hand size
 * vector. The permutation information are required for the solving.
 *
 * L and U matrix can safely alias each other since they are read only. In
 * fact, their subelements do not alias each other but only stored in the
 * same matrix A.
 */
static int FactorizedLinearSystemSolver(const int n, Real L[restrict][n], Real U[restrict][n], 
        Real x[], Real b[], const int permute[restrict])
{
    /*
     * Rearrange the elements of right hand side vector according to
     * permutations applied. Rearranged vector can be saved to the 
     * solution vector, interchange values needs to be applied since 
     * the right hand side vector and solution vector may occupy the
     * same storage space.
     */
    Real temp = 0.0; /* auxiliary variable */
    for (int row = 0; row < n; ++row) {
        temp = b[row];
        x[row] = b[permute[row]];
        b[permute[row]] = temp;
    }
    /*
     * Forward substitution.
     * Note that the lower triangle L has unit diagonal.
     */
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < row; ++col) {
            x[row] = x[row] - L[row][col] * x[col];
        }
    }
    /*
     * Bach substitution.
     */
    for (int row = n - 1; row >= 0; --row) {
        for (int col = row + 1; col < n; ++col) {
            x[row] = x[row] - U[row][col] * x[col];
        }
        x[row] = x[row] / U[row][row];
    }
    return 0;
}
/* a good practice: end file with a newline */

