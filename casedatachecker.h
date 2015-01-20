/****************************************************************************
 * Header File                                                              *
 * Last-modified: 18 Jan 2015 09:31:56 PM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_CASEDATACHECKER_H_ /* if this is the first definition */
#define ARTRACFD_CASEDATACHECKER_H_ /* a unique marker for this header file */
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
int CheckCaseSettingData(const Space *, const Time *, const Fluid *,
        const Reference *);
#endif
/* a good practice: end file with a newline */

