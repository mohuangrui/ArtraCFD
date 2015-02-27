/****************************************************************************
 * Numeric Scheme for Space Domain                                          *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines the numeric schemes of space domain.                 *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "tvd.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static Real Q(const Real x);
static int sgn(const Real x);
static Real minmod(const Real x, const Real y, const Real z);
static Real Min(const Real valueA, const Real valueB);
static Real Max(const Real valueA, const Real valueB);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int TVD(Real *fieldDataU, const Real *fieldDataUn, const Space *space, 
        const Partition *part, const Flow *flow)
{
    /*
     * Decompose the conservative field variable into each component.
     */
    Real *U[5] = {
        fieldDataU + 0 * space->nMax,
        fieldDataU + 1 * space->nMax,
        fieldDataU + 2 * space->nMax,
        fieldDataU + 3 * space->nMax,
        fieldDataU + 4 * space->nMax};
    const Real *Un[5] = {
        fieldDataUn + 0 * space->nMax,
        fieldDataUn + 1 * space->nMax,
        fieldDataUn + 2 * space->nMax,
        fieldDataUn + 3 * space->nMax,
        fieldDataUn + 4 * space->nMax};
    /*
     * Decompose the nonviscous flux variables into each component
     */
    Real Fx[5] = {0, 0, 0, 0, 0};
    Real Fy[5] = {0, 0, 0, 0, 0};
    Real Fz[5] = {0, 0, 0, 0, 0};
    /*
     * Decompose the viscous flux variables into each component
     */
    Real Gx[5] = {0, 0, 0, 0, 0};
    Real Gy[5] = {0, 0, 0, 0, 0};
    Real Gz[5] = {0, 0, 0, 0, 0};
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
    int idxW = 0; /* index at West */
    int idxE = 0; /* index at East */
    int idxS = 0; /* index at South */
    int idxN = 0; /* index at North */
    int idxF = 0; /* index at Front */
    int idxB = 0; /* index at Back */
    /*
     * When exchange a large bunch of data between two arrays, if there is no
     * new data generation but just data exchange and update, then the rational
     * way is to exchange the head address that they  points rather than values
     * of data entries.
     */
    return 0;
}
static int Lx(Real *fieldDataU, const Real *fieldDataUn, const Space *space, 
        const Partition *part, const Flow *flow)
{
    /*
     * Decompose the conservative field variable into each component.
     */
    Real *U[5] = {
        fieldDataU + 0 * space->nMax,
        fieldDataU + 1 * space->nMax,
        fieldDataU + 2 * space->nMax,
        fieldDataU + 3 * space->nMax,
        fieldDataU + 4 * space->nMax};
    const Real *Un[5] = {
        fieldDataUn + 0 * space->nMax,
        fieldDataUn + 1 * space->nMax,
        fieldDataUn + 2 * space->nMax,
        fieldDataUn + 3 * space->nMax,
        fieldDataUn + 4 * space->nMax};
    /*
     * Decompose the nonviscous flux variables into each component
     */
    Real Fx[5] = {0, 0, 0, 0, 0};
    /*
     * Decompose the viscous flux variables into each component
     */
    Real Gx[5] = {0, 0, 0, 0, 0};
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
    int idxW = 0; /* index at West */
    int idxE = 0; /* index at East */
    return 0;
}
static Real Q(const Real x)
{
    const Real e = 0.01;
    if (fabs(x) >= e) {
        return fabs(x);
    }
    return (0.5 * (x * x / e + e));
}
static int sgn(const Real x)
{
    if (x > 0) {
        return 1;
    }
    if (x < 0) {
        return -1;
    }
    return 0;
}
static Real minmod(const Real x, const Real y, const Real z)
{
    if ((x * y <= 0) || (x * z <= 0)) {
        return 0;
    }
    return (sgn(x) * Min(fabs(x), Min(fabs(y), fabs(z))));
}
static Real Min(const Real valueA, const Real valueB)
{
    if (valueA < valueB) {
        return valueA;
    }
    return valueB;
}
static Real Max(const Real valueA, const Real valueB)
{
    if (valueA > valueB) {
        return valueA;
    }
    return valueB;
}
/* a good practice: end file with a newline */

