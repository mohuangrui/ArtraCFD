/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_POSTPROCESS_H_ /* if this is the first definition */
#define ARTRACFD_POSTPROCESS_H_ /* a unique marker for this header file */
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
 * Postprocessor
 *
 * Function
 *      Perform postprocessing for ArtraCFD program, such as memory release.
 */
int Postprocess(Field *, Flux *, Space *, Particle *, Partition *);
#endif
/* a good practice: end file with a newline */

