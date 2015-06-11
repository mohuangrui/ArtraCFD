/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_DOMAINPARTITION_H_ /* if this is the first definition */
#define ARTRACFD_DOMAINPARTITION_H_ /* a unique marker for this header file */
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
 * Domain partition
 *
 * Function
 *      Decomposite the computational domain, and provide the index range
 *      for each partitions, as well as partition names.
 */
extern int DomainPartition(Partition *, const Space *);
#endif
/* a good practice: end file with a newline */

