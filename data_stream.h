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
#ifndef ARTRACFD_DATA_STREAM_H_ /* if undefined */
#define ARTRACFD_DATA_STREAM_H_ /* set a unique marker */
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
extern void WriteData(const int n, const Time *, const Space *, const Model *);
extern void ReadData(const int n, Time *, Space *, const Model *);
extern void WritePolyStateData(const int pm, const int pn, FILE *fp, const Geometry *const);
extern void ReadPolyStateData(const int pm, const int pn, FILE *fp, Geometry *const);
#endif
/* a good practice: end file with a newline */

