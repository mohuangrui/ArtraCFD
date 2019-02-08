/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>                        *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_DATA_PROBE_H_ /* if undefined */
#define ARTRACFD_DATA_PROBE_H_ /* set a unique marker */
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
extern void WritePointProbeData(const Time *, const Space *, const Model *);
extern void WriteLineProbeData(const Time *, const Space *, const Model *);
extern void WriteCurveProbeData(const Time *, const Space *, const Model *);
extern void WriteSurfaceForceData(const Time *, const Space *, const Model *);
#endif
/* a good practice: end file with a newline */

