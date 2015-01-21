/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_PREPROCESS_H_ /* if this is the first definition */
#define ARTRACFD_PREPROCESS_H_ /* a unique marker for this header file */
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
 * Preprocessor
 *
 * Parameter
 *      data structures
 * Function
 *      Call a series of functions to perform the preprocessing of
 *      computatioal fluid dynamics.
 * Returns
 *      0 -- successful
 */
int Preprocess(Field *, Flux *, Space *, Particle *, Time *,
        Partition *, Fluid *, Flow *, Reference *);
#endif
/* a good practice: end file with a newline */

