/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_GEOMETRYLOADER_H_ /* if this is the first definition */
#define ARTRACFD_GEOMETRYLOADER_H_ /* a unique marker for this header file */
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
 * Geometry data loader
 *
 * Function
 *      Load geometry data to program.
 *      -- non restart run, load from file artracfd.geo
 *      -- restart run, load from the restart particle file.
 */
int LoadGeometryData(Particle *, const Time *);
#endif
/* a good practice: end file with a newline */

