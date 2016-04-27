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
/*
 * Function pointers are useful for implementing a form of polymorphism.
 * They are mainly used to reduce or avoid switch statement. Pointers to
 * functions can get rather messy. Declaring a typedef to a function pointer
 * generally clarifies the code.
 */
typedef int (*StructuredDataWriter)(const Space *, const Time *, const Model *);
typedef int (*StructuredDataReader)(Space *, Time *, const Model *);
typedef int (*PolyDataWriter)(const Geometry *, const Time *);
typedef int (*PolyDataReader)(Geometry *, const Time *);
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
int WriteFieldData(const Space *space, const Time *time, const Model *model)
{
    WriteStructuredData[time->dataStreamer](space, time, model);
    return 0;
}
int ReadFieldData(Space *space, Time *time, const Model *model)
{
    ReadStructuredData[time->dataStreamer](space, time, model);
    return 0;
}
int WriteGeometryData(const Geometry * geo, const Time *time)
{
    if (0 != geo->totalN) {
        WritePolyData[0](geo, time);
    }
    return 0;
}
int ReadGeometryData(Geometry * geo, const Time *time)
{
    if (0 != geo->totalN) {
        ReadPolyData[0](geo, time);
    }
    return 0;
}
/* a good practice: end file with a newline */

