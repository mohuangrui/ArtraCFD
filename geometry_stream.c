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
#include "geometry_stream.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "stl.h"
#include "paraview.h"
#include "computational_geometry.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NonrestartGeometryReader(Geometry *);
static int RestartGeometryReader(const Time *, Geometry *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ReadGeometryData(Space *space, const Time *time)
{
    if (0 == time->restart) {
        NonrestartGeometryReader(&(space->geo));
    } else {
        RestartGeometryReader(time, &(space->geo));
    }
    ComputeGeometryParameters(space, &(space->geo));
    return 0;
}
static int NonrestartGeometryReader(Geometry *geo)
{
    ShowInformation("Reading geometry data ...");
    FILE *filePointer = fopen("artracfd.geo", "r");
    if (NULL == filePointer) {
        FatalError("failed to open file: artracfd.geo...");
    }
    /* read and process file line by line */
    String currentLine = {'\0'}; /* store the current read line */
    String fileName = {'\0'}; /* store current geometry file name */
    int entryCount = 0; /* entry count */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "count begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(geo->totalN));
            if (0 >= geo->totalN) {
                geo->totalN = 0;
                ++entryCount;
                break;
            }
            geo->list = AssignStorage(geo->totalN, "Polyhedron");
            continue;
        }
        if (0 == strncmp(currentLine, "sphere begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(geo->sphereN));
            if (0 >= geo->sphereN) {
                geo->sphereN = 0;
                continue;
            }
            if (geo->totalN < geo->sphereN) {
                geo->sphereN = geo->totalN;
            }
            fgets(currentLine, sizeof currentLine, filePointer);
            
            for (int n = 0; n < geo->sphereN; ++n) {
                ReadPolyhedronStateData(&filePointer, geo->list + n);
                geo->list[n].facetN = 0; /* analytical geometry tag */
                geo->list[n].facet = NULL;
            }
            if (geo->totalN == geo->sphereN) {
                break;
            }
            continue;
        }
        if (0 == strncmp(currentLine, "STL begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%s", fileName);
            ReadStlFile(fileName, geo->list + geo->sphereN + geo->stlN);
            ReadPolyhedronStateData(&filePointer, geo->list + geo->sphereN + geo->stlN);
            ++geo->stlN; /* point to the next geometry */
            continue;
        }
    }
    fclose(filePointer); /* close current opened file */
    /* Check missing information section in configuration */
    if (2 + geo->stlN != entryCount) {
        FatalError("missing necessary information section");
    }
    ShowInformation("Session End");
    return 0;
}
static int RestartGeometryReader(const Time *time, Geometry *geo)
{
    ShowInformation("Restore geometry data ...");
    ReadGeometryDataParaview(time, geo);
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

