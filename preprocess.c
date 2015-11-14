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
#include "geometry_stream.h"
#include "immersed_boundary_treatment.h"
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
int Preprocess(Space *space, Time *time, Model *model, Partition *part, Geometry *geo)
{
    LoadCaseSettingData(space, time, model, part);
    ComputeCFDParameters(space, time, model);
    DomainPartition(space, part);
    ProgramMemoryAllocate(space);
    ReadGeometryData(space, time, model, geo);
    ComputeGeometryDomain(space, part, geo);
    return 0;
}
/*
 * This function allocates memory for field data. The storage retrieving
 * need to be done in the postprocessor.
 */
static int ProgramMemoryAllocate(Space *space)
{
    ShowInformation("Allocating memory...");
    space->node = AssignStorage(space->totalN, "Node");
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

