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
 * Function definitions
 ****************************************************************************/
int WriteComputedData(const Real *U, const Space *space, const Time *time,
        const Model *model, const Partition *part)
{
    switch (time->dataStreamer) {
        case 0: /* ParaView */
            WriteComputedDataParaview(U, space, time, model, part);
            break;
        case 1: /* Generic Ensight */
            WriteComputedDataEnsight(U, space, time, model, part);
            break;
        case 2: /* ParaView Ensight */
            WriteComputedDataParasight(U, space, time, model, part);
            break;
        default: /* ParaView */
            WriteComputedDataParaview(U, space, time, model, part);
            break;
    }
    return 0;
}
int LoadComputedData(Real *U, const Space *space, Time *time,
        const Model *model, const Partition *part)
{
    switch (time->dataStreamer) {
        case 0: /* ParaView */
            LoadComputedDataParaview(U, space, time, model, part);
            break;
        case 1: /* Ensight */
            LoadComputedDataEnsight(U, space, time, model, part);
            break;
        case 2: /* Paraview Ensight */
            LoadComputedDataParasight(U, space, time, model, part);
            break;
        default: /* ParaView */
            LoadComputedDataParaview(U, space, time, model, part);
            break;
    }
    return 0;
}
/* a good practice: end file with a newline */

 
