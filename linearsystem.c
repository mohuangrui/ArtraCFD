/****************************************************************************
 * Linear System Solver                                                     *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int LUFactorization(Real **A, const int n, int permute[]);
static int FactorizedLinearSystemSolver(Real **L, Real **U, Real x[], Real b[],
        const int n, const int permute[]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int MatrixLinearSystemSolver(Real **A, Real **X, Real **B, const int n)
{
    int permute[n]; /* record the permutation information */
    Real rhs[n]; /* transfer data into column vector */
    LUFactorization(A, n, permute);
    for (int col = 0; col < n; ++col) { /* solve column by column */
        for (int row = 0; row < n; ++row) {
            rhs[row] = B[row][col]; /* obtain each right hand side vector */
        }
        FactorizedLinearSystemSolver(A, A, rhs, rhs, n, permute);
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
 */
static int LUFactorization(Real **A, const int n, int permute[])
{
    const Real epsilon = 1.0e-15; /* a small number for singularity check */
    Real temp = 0.0; /* auxiliary variable */
    Real scale[n]; /* stores the implicit scaling of each row with variable length array */
    int rowMax = 0; /* record the row number of current pivot element */
    int sign = 1; /* used for determine the sign of determinant after pivoting */
    /*
     * Loop over rows to get the implicit scaling information.
     */
    for (int row = 0; row < n; ++row) {
        Real maximum = 0.0;
        for (int col = 0; col < n; ++col) {
            temp = fabs(A[row][col]);
            if (temp > maximum) {
                maximum = temp;
            }
        }
        if (0 == maximum) {
            FatalError("singular matrix in LU factorization...");
        }
        scale[row] = 1.0 / maximum; /* save the scaling */
    }
    /*
     * Do LU factorization with partial pivoting.
     */
    for (int loop = 0; loop < n; ++loop) {
        Real maximum = 0.0; /* initialize for the search of largest pivot element */
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
        if (0 == A[loop][loop]) { /* substitute zero pivot element with epsilon */
            A[loop][loop] = epsilon;
        }
        for (int row = loop + 1; row < n; ++row) {
            A[row][loop] = A[row][loop] / A[loop][loop]; /* divide column element with pivot element */
            for (int col = loop + 1; col < n; ++col) {
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
 */
static int FactorizedLinearSystemSolver(Real **L, Real **U, Real x[], Real b[],
        const int n, const int permute[])
{
    /*
     * Rearrange the elements of right hand side vector according to
     * permutations applied. Rearranged vector can be saved to the 
     * solution vector, interchange values needs to be applied since 
     * the right hand side vector and solution vector may occupy the
     * same storage space.
     */
    for (int row = 0; row < n; ++row) {
        const Real temp = b[row];
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
}
/* a good practice: end file with a newline */

