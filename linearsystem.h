/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_LINEARSYSTEM_H_ /* if this is the first definition */
#define ARTRACFD_LINEARSYSTEM_H_ /* a unique marker for this header file */
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
 *      A is the system matrix, X is the solution matrix with each column vector
 *      stands for a solution vector of a single linear system A * x = b, B is 
 *      the right hand matrix with each column vector stands for a right hand 
 *      vector of a single linear system. If B is a unit matrix, then solution
 *      matrix X is the inverse matrix of A. In practice, right hand matrix 
 *      and solution matrix can share the same storage.
 */
extern int MatrixLinearSystemSolver(Real **A, Real **X, Real **B, const int n);
#endif
/* a good practice: end file with a newline */

