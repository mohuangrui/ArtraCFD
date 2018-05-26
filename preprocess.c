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
    ShowInformation("Session End");
    ShowInformation("Preprocessing...");
    fprintf(stdout, "  loading case data...\n");
    LoadCaseData(time, space, model);
    fprintf(stdout, "  computing parameters...\n");
    ComputeParameters(time, space, model);
    fprintf(stdout, "  partitioning domain...\n");
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
    Partition *part = &(space->part);
    Geometry *geo = &(space->geo);
    const int totN = part->n[X] * part->n[Y] * part->n[Z];
    space->node = AssignStorage(totN * sizeof(*space->node));
    if (0 != geo->totN) {
        geo->col = AssignStorage(geo->totN * sizeof(*geo->col));
        geo->poly = AssignStorage(geo->totN * sizeof(*geo->poly));
    }
    return 0;
}
/* a good practice: end file with a newline */

