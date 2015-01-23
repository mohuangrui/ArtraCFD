/****************************************************************************
 * Header File                                                              *
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
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Ensight transient case file initializer
 *
 * Function
 *      Initialize a Ensight transient case file. This function only needs to
 *      be call once for each non restart run.
 */
int InitializeEnsightTransientCaseFile(const Time *);
/*
 * Ensight format data exporter
 *
 * Function
 *      Export primitive field data vector variable
 *      fieldData = [rho, u, v, w, p, T]
 *      to binary data files with Ensight data format.
 * Notice
 *      fieldData is a linear array that stores all the values of all the 
 *      primitive field variables. These data are in sequential state and can
 *      be accessed by linear index math.
 */
int WriteComputedDataEnsight(const Real * fieldData, const Space *, 
        const Particle *, const Time *, const Partition *);
#endif
/* a good practice: end file with a newline */

