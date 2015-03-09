/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_TEMPERALMARCHING_H_ /* if this is the first definition */
#define ARTRACFD_TEMPERALMARCHING_H_ /* a unique marker for this header file */
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
 * Temperal marching
 *
 * Function
 */
extern int TemperalMarching(Field *, Space *, Particle *, Time *, 
        const Partition *, const Flow *);
#endif
/* a good practice: end file with a newline */

