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
#include "postprocess.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ProgramMemoryRelease(Space *);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
int Postprocess(Space *space)
{
    ShowInformation("Postprocessing...");
    fprintf(stdout, "  releasing memory...\n");
    ProgramMemoryRelease(space);
    fprintf(stdout, "  computing finished, successfully exit!n");
    ShowInformation("Session End");
    return 0;
}
static int ProgramMemoryRelease(Space *space)
{
    /* geometry related */
    for (int n = space->geo.sphereN; n < space->geo.totalN; ++n) {
        RetrieveStorage(space->geo.list[n].facet);
    }
    RetrieveStorage(space->geo.list);
    /* field variable related */
    RetrieveStorage(space->node);
    return 0;
}
/* a good practice: end file with a newline */

