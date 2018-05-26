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
 * Required Header Files
 ****************************************************************************/
#include "data_stream.h"
#include <stdio.h> /* standard library for input and output */
#include "paraview.h"
#include "ensight.h"
#include "commons.h"
/****************************************************************************
 * Function Pointers
 ****************************************************************************/
typedef int (*StructuredDataWriter)(const Time *, const Space *, const Model *);
typedef int (*StructuredDataReader)(Time *, Space *, const Model *);
typedef int (*PolyDataWriter)(const Time *, const Geometry *);
typedef int (*PolyDataReader)(const Time *, Geometry *);
/****************************************************************************
 * Global Variables Definition with Private Scope
 ****************************************************************************/
static StructuredDataWriter WriteStructuredData[2] = {
    WriteStructuredDataParaview,
    WriteStructuredDataEnsight};
static StructuredDataReader ReadStructuredData[2] = {
    ReadStructuredDataParaview,
    ReadStructuredDataEnsight};
static PolyDataWriter WritePolyData[1] = {
    WritePolyDataParaview};
static PolyDataReader ReadPolyData[1] = {
    ReadPolyDataParaview};
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteFieldData(const Time *time, const Space *space, const Model *model)
{
    WriteStructuredData[time->dataStreamer](time, space, model);
    return 0;
}
int ReadFieldData(Time *time, Space *space, const Model *model)
{
    ReadStructuredData[time->dataStreamer](time, space, model);
    return 0;
}
int WriteGeometryData(const Time *time, const Geometry *geo)
{
    if (0 == geo->totN) {
        return 0;
    }
    WritePolyData[0](time, geo);
    return 0;
}
int ReadGeometryData(const Time *time, Geometry *geo)
{
    if (0 == geo->totN) {
        return 0;
    }
    ReadPolyData[0](time, geo);
    return 0;
}
/* a good practice: end file with a newline */

