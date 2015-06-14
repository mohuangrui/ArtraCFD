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
 * Required Header Files
 ****************************************************************************/
#include "data_stream.h"
#include <stdio.h> /* standard library for input and output */
#include "paraview_stream.h"
#include "parasight_stream.h"
#include "ensight_stream.h"
#include "commons.h"
/****************************************************************************
 * Function pointers
 ****************************************************************************/
/*
 * Function pointers are useful for implementing a form of polymorphism.
 * They are mainly used to reduce or avoid switch statement. Pointers to
 * functions can get rather messy. Declaring a typedel to a function pointer
 * generally clarifies the code.
 */
typedef int (*DataWriter)(const Real *, const Space *, const Time *,
        const Model *, const Partition *);
typedef int (*DataReader)(Real *, const Space *, Time *,
        const Model *, const Partition *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteComputedData(const Real *U, const Space *space, const Time *time,
        const Model *model, const Partition *part)
{
    DataWriter WriteData[3] = {
        WriteComputedDataParaview,
        WriteComputedDataEnsight,
        WriteComputedDataParasight
    };
    WriteData[time->dataStreamer](U, space, time, model, part);
    return 0;
}
int ReadComputedData(Real *U, const Space *space, Time *time,
        const Model *model, const Partition *part)
{
    DataReader ReadData[3] = {
        ReadComputedDataParaview,
        ReadComputedDataEnsight,
        ReadComputedDataParasight
    };
    ReadData[time->dataStreamer](U, space, time, model, part);
    return 0;
}
/* a good practice: end file with a newline */

