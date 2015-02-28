/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_INITIALIZATION_H_ /* if this is the first definition */
#define ARTRACFD_INITIALIZATION_H_ /* a unique marker for this header file */
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
 * Flow field initializer
 *
 * Function
 *      Initialize the flow field variable.
 *      non restart -- initialize by specified data.
 *      restart -- initialize by restart data.
 */
extern int InitializeFlowField(Real *U, const Space *, const Particle *, Time *,
        const Partition *, const Flow *);
#endif
/* a good practice: end file with a newline */

