/****************************************************************************
 * Header File                                                              *
 * Last-modified: 17 Jan 2015 12:22:47 PM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_SOLVE_H_ /* if this is the first definition */
#define ARTRACFD_SOLVE_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/****************************************************************************
 * Function declaration
 ****************************************************************************/
int Solve(Field *, Flux *, Space *, Particle *, Time *, 
        const Partition *, const Fluid *, const Flow *, const Reference *);
#endif
/* a good practice: end file with a newline */

