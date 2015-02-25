/****************************************************************************
 * Header File                                                              *
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
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Initialize domain geometry
 *
 * Function
 *      Initialize domain geometry. It only needs to be call once.
 */
extern int InitializeDomainGeometry(Space *);
/*
 * Compute domain geometry
 *
 * Function
 *      Employ GCIBM approach to handle complex geometry that locates in
 *      the computational domain.
 */
extern int ComputeDomainGeometryGCIBM(Space *, Particle *, const Partition *);
#endif
/* a good practice: end file with a newline */

