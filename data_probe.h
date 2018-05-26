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
#ifndef ARTRACFD_DATA_PROBE_H_ /* if this is the first definition */
#define ARTRACFD_DATA_PROBE_H_ /* a unique marker for this header file */
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
 * Write probe data
 *
 * Function
 *      Write field data of probes.
 */
extern void WriteFieldDataAtPointProbes(const Time *, const Space *, const Model *);
extern void WriteFieldDataAtLineProbes(const Time *, const Space *, const Model *);
extern void WriteFieldDataAtCurveProbes(const Time *, const Space *, const Model *);
extern void WriteSurfaceForceData(const Time *, const Space *);
#endif
/* a good practice: end file with a newline */

