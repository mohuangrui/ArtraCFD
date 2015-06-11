/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) 2014-2018 Huangrui Mo <huangrui.mo@gmail.com>              *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_GEOMETRY_STREAM_H_ /* if this is the first definition */
#define ARTRACFD_GEOMETRY_STREAM_H_ /* a unique marker for this header file */
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
 *      -- restart run, load from the restart geometry file.
 */
extern int LoadGeometryData(Geometry *, const Space *, const Time *, const Flow *);
/*
 * Geometry data writer
 *
 * Function
 *      Write geometry data.
 */
extern int WriteGeometryData(const Geometry *, const Time *);
#endif
/* a good practice: end file with a newline */

