/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_BOUNDARYCONDITION_H_ /* if this is the first definition */
#define ARTRACFD_BOUNDARYCONDITION_H_ /* a unique marker for this header file */
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
 * Boundary condition
 *
 * Function
 *      Apply boundary conditions and treatments for the flow field variable.
 */
extern int BoundaryCondtionsAndTreatments(Real *U, const Space *, const Particle *,
        const Partition *, const Flow *);
#endif
/* a good practice: end file with a newline */

