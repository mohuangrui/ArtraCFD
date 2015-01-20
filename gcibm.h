/****************************************************************************
 * Header File                                                              *
 * Last-modified: 18 Jan 2015 10:06:14 PM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_GCIBM_H_ /* if this is the first definition */
#define ARTRACFD_GCIBM_H_ /* a unique marker for this header file */
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
int InitializeDomainGeometryGCIBM(Space *, Particle *, const Partition *);
int ComputeDomainGeometryGCIBM(Space *, Particle *, const Partition *);
#endif
/* a good practice: end file with a newline */

