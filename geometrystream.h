/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_GEOMETRYSTREAM_H_ /* if this is the first definition */
#define ARTRACFD_GEOMETRYSTREAM_H_ /* a unique marker for this header file */
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
extern int LoadGeometryData(Particle *, const Time *);
/*
 * Geometry data writer
 *
 * Function
 *      Write geometry data.
 */
extern int WriteGeometryData(const Particle *, const Time *);
#endif
/* a good practice: end file with a newline */

