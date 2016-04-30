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
#include "preprocess.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include "case_loader.h"
#include "cfd_parameters.h"
#include "domain_partition.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ProgramMemoryAllocate(Space *);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * This is the overall preprocessing function
 */
int Preprocess(Time *time, Space *space, Model *model)
{
    ShowInformation("Preprocessing...");
    fprintf(stdout, "  loading case setting data...\n");
    LoadCaseSettingData(time, space, model);
    fprintf(stdout, "  computing parameters...\n");
    ComputeCFDParameters(time, space, model);
    fprintf(stdout, "  domain partitioning...\n");
    DomainPartition(space);
    fprintf(stdout, "  allocating memory...\n");
    ProgramMemoryAllocate(space);
    ShowInformation("Session End");
    return 0;
}
/*
 * This function allocates memory for field data. The storage retrieving
 * need to be done in the postprocessor.
 */
static int ProgramMemoryAllocate(Space *space)
{
    space->node = AssignStorage(space->part.n[X] * space->part.n[Y] * space->part.n[Z], "Node");
    if (0 != space->geo.totalN) {
        space->geo.poly = AssignStorage(space->geo.totalN, "Polyhedron");
    }
    return 0;
}
/* a good practice: end file with a newline */

