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
#endif
/* a good practice: end file with a newline */

