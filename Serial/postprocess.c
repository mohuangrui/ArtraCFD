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
static int ProgramMemoryRelease(Time *, Space *);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
int Postprocess(Time *time, Space *space)
{
    ShowInformation("Postprocessing...");
    fprintf(stdout, "  releasing memory...\n");
    ProgramMemoryRelease(time, space);
    fprintf(stdout, "  computing finished, successfully exit.\n");
    ShowInformation("Session End");
    return 0;
}
static int ProgramMemoryRelease(Time *time, Space *space)
{
    /* geometry related */
    Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    for (int n = geo->sphN; n < geo->totN; ++n) {
        poly = geo->poly + n;
        RetrieveStorage(poly->f);
        RetrieveStorage(poly->Nf);
        RetrieveStorage(poly->e);
        RetrieveStorage(poly->Ne);
        RetrieveStorage(poly->v);
        RetrieveStorage(poly->Nv);
    }
    RetrieveStorage(geo->poly);
    RetrieveStorage(geo->col);
    /* field variable related */
    RetrieveStorage(space->node);
    /* time related */
    RetrieveStorage(time->lp);
    RetrieveStorage(time->pp);
    return 0;
}
/* a good practice: end file with a newline */

