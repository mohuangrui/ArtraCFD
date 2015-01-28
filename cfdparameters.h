/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_CFDPARAMETERS_H_ /* if this is the first definition */
#define ARTRACFD_CFDPARAMETERS_H_ /* a unique marker for this header file */
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
 * Compute CFD parameters
 *
 * Function
 *      Compute parameters required for CFD, and perform nondimensionalize 
 *      operations on fluid and flow variables to unify the dimensional 
 *      form and nondimensional form of governing equations.
 */
extern int ComputeCFDParameters(Space *, Time *, const Fluid *, Flow *,
        const Reference *);
#endif
/* a good practice: end file with a newline */

