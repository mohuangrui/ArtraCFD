/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_CFDCOMMONS_H_ /* if this is the first definition */
#define ARTRACFD_CFDCOMMONS_H_ /* a unique marker for this header file */
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
 * Compute nonviscous flux variables
 *
 * Function
 *      Compute flux variables based on conservative variables.
 */
extern int ComputeNonviscousFlux(const Field *, Flux *, const Space *, const Flow *);
/*
 * Compute viscous flux variables
 *
 * Function
 *      Compute flux variables based on conservative variables.
 */
extern int ComputeViscousFlux(const Field *, Flux *, const Space *, const Flow *);
/*
 * Compute dt by CFL condition
 *
 * Function
 *      Compute the time step size dt based on CFL condition.
 *
 * Returns
 *      a real -- the calculated time step size dt.
 */
extern Real ComputeTimeStepByCFL(const Field *, const Space *, const Time *, 
        const Partition *, const Flow *);
/*
 * Minimum positive
 *
 * Function
 *      Find the minimum positive value among two values, if none of them is
 *      positive, return a very large number instead.
 *
 * Returns
 *      a real -- the smaller value which is positive
 *      1e38 -- if none of them is positive
 */
extern Real MinPositive(const Real valueA, const Real valueB);
/*
 * Minimum
 *
 * Function
 *      Find the minimum value among two values.
 *
 * Returns
 *      a real -- the smaller value
 */
extern Real Min(const Real valueA, const Real valueB);
/*
 * Maximum
 *
 * Function
 *      Find the maximum value among two values.
 *
 * Returns
 *      a real -- the larger value
 */
extern Real Max(const Real valueA, const Real valueB);
#endif
/* a good practice: end file with a newline */

