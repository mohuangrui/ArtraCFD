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
static void AllocateProgramMemory(Space *, Model *);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
int Preprocess(Time *time, Space *space, Model *model)
{
    ShowInfo("Session");
    ShowInfo("Preprocessing...\n");
    ShowInfo("  loading case data...\n");
    LoadCaseData(time, space, model);
    ShowInfo("  computing parameters...\n");
    ComputeParameters(time, space, model);
    ShowInfo("  partitioning domain...\n");
    PartitionDomain(space);
    ShowInfo("  allocating memory...\n");
    AllocateProgramMemory(space, model);
    ShowInfo("Session");
    return 0;
}
/*
 * Allocate memory for the remaining unassigned data.
 * Storage retrieving is done in the postprocessor.
 */
static void AllocateProgramMemory(Space *space, Model *model)
{
    Partition *const part = &(space->part);
    Geometry *const geo = &(space->geo);
    const int totN = part->n[X] * part->n[Y] * part->n[Z];
    space->node = AssignStorage(totN * sizeof(*space->node));
    if (0 != geo->totN) {
        geo->col = AssignStorage(geo->totN * sizeof(*geo->col));
        geo->poly = AssignStorage(geo->totN * sizeof(*geo->poly));
    }
    model->mat = AssignStorage(sizeof(*model->mat));
    return;
}
/* a good practice: end file with a newline */

