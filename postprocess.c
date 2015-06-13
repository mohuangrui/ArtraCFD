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
#include "postprocess.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ProgramMemoryRelease(Field *, Space *, Geometry *);
static int FinalInformation(void);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * This is the overall postprocessing function
 */
int Postprocess(Field *field, Space *space, Geometry *geometry)
{
    ProgramMemoryRelease(field, space, geometry);
    FinalInformation();
    return 0;
}
/*
 * This function together with some subfuctions realize the dynamic
 * memory release for each global pointer
 */
static int ProgramMemoryRelease(Field *field, Space *space, Geometry *geometry)
{
    ShowInformation("Releasing memory back to system...");
    /* field variable related */
    RetrieveStorage(field->Un);
    RetrieveStorage(field->U);
    RetrieveStorage(field->Uswap);
    /* space related */
    RetrieveStorage(space->nodeFlag);
    /* geometry related */
    RetrieveStorage(geometry->headAddress);
    ShowInformation("Session End");
    return 0;
}
/*
 * This function shows some final information
 */
static int FinalInformation(void)
{
    ShowInformation("Computing finished, successfully exit!");
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

 
