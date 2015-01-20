/****************************************************************************
 * Header File                                                              *
 * Last-modified: 17 Jan 2015 12:38:13 PM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_ENSIGHT_H_ /* if this is the first definition */
#define ARTRACFD_ENSIGHT_H_ /* a unique marker for this header file */
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
int InitializeEnsightTransientCaseFile(const Time *);
int WriteComputedDataEnsight(const double *, const Space *, const Particle *,
        const Time *, const Partition *);
#endif
/* a good practice: end file with a newline */

