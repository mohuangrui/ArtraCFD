/****************************************************************************
 * Header File                                                              *
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
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Case setting data loader
 *
 * Function
 *      Load case setting data from file artracfd.case, if file does not 
 *      exist or any illegal data exists in the file, it terminates the
 *      program.
 * Outcome
 *      Will write loaded data to file artracfd.verify for verification.
 */
extern int LoadCaseSettingData(Space *, Time *, Fluid *, Reference *);
#endif
/* a good practice: end file with a newline */

