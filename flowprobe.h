/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_FLOWPROBE_H_ /* if this is the first definition */
#define ARTRACFD_FLOWPROBE_H_ /* a unique marker for this header file */
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
 * Write probe data
 *
 * Function
 *      Write field data of probes.
 */
extern int WriteComputedDataAtProbes(const int stepCount, const Real *U, 
        const Space *, const Flow *);
#endif
/* a good practice: end file with a newline */

