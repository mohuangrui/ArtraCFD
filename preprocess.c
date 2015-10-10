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
#include "preprocess.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include "case_loader.h"
#include "cfd_parameters.h"
#include "domain_partition.h"
#include "geometry_stream.h"
#include "gcibm.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ProgramMemoryAllocate(Field *, Space *);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * This is the overall preprocessing function
 */
int Preprocess(Field *field, Space *space, Time *time, Model *model,
        Partition *part, Geometry *geometry)
{
    LoadCaseSettingData(space, time, model, part);
    ComputeCFDParameters(space, time, model);
    DomainPartition(space, part);
    ProgramMemoryAllocate(field, space);
    LoadGeometryData(space, time, model, geometry);
    ComputeGeometryDomain(space, part, geometry);
    return 0;
}
/*
 * This function together with some subfuctions realize the dynamic
 * memory allocation for each global pointer. The storage retrieving
 * need to be done in the postprocessor.
 */
static int ProgramMemoryAllocate(Field *field, Space *space)
{
    ShowInformation("Allocating memory...");
    /*
     * Conservative variables.
     * Tips: the storage space of U is best between Un and Uswap.
     */
    int idxMax = space->nMax * DIMU;
    field->Un = AssignStorage(idxMax, "Real");
    field->U = AssignStorage(idxMax, "Real");
    field->Uswap = AssignStorage(idxMax, "Real");
    /*
     * Node type identifier
     */
    idxMax = space->nMax;
    space->nodeFlag = AssignStorage(idxMax, "int");
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

