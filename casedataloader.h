/****************************************************************************
 * Header File                                                              *
 * Last-modified: 18 Jan 2015 09:16:13 PM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_CASEDATALOADER_H_ /* if this is the first definition */
#define ARTRACFD_CASEDATALOADER_H_ /* a unique marker for this header file */
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
int LoadCaseSettingData(Space *, Time *, Fluid *, Reference *);
#endif
/* a good practice: end file with a newline */

