/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_RUNGEKUTTA_H_ /* if this is the first definition */
#define ARTRACFD_RUNGEKUTTA_H_ /* a unique marker for this header file */
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
 * Time marching scheme
 *
 * Function
 *      Apply runge-kutta method for explicit time marching solution
 */
extern int RungeKuttaTimeMarching(Field *, Space *, Particle *, Time *, 
        const Partition *, const Flow *);
#endif
/* a good practice: end file with a newline */

