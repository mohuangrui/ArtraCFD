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
#include "case_loader.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ReadCaseSettingData(Space *, Time *, Model *, Partition *);
static int ReadBoundaryData(FILE **, Partition *, const int);
static int ReadConsecutiveRealData(FILE **, Real *, const int);
static int WriteBoundaryData(FILE **, const Partition *, const int);
static int WriteRegionalInitializerData(FILE **, const Partition *, const int);
static int WriteVerifyData(const Space *, const Time *, const Model *, const Partition *);
static int CheckCaseSettingData(const Space *, const Time *, const Model *, const Partition *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function load the case settings from the case file.
 */
int LoadCaseSettingData(Space *space, Time *time, Model *model, Partition *part)
{
    ShowInformation("Loading case setting data ...");
    ReadCaseSettingData(space, time, model, part);
    WriteVerifyData(space, time, model, part);
    CheckCaseSettingData(space, time, model, part);
    ShowInformation("Session End");
    return 0;
}
/*
 * This function read the case settings from the case file.
 * The key is to read and process file line by line. Use "*** begin" 
 * in the case file to identify and control the reading. 
 * The function scanf is notorious for its poor end-of-line handling.
 * Instead, use fgets to read a line of input and sscanf to process it.
 * Note: use a large enough number when using fgets to ensure reading
 * a whole line at a time. fgets will get the entire line including
 * the newline character (\n).
 * NOTE: if memory locations of input objects overlap, the behavior of
 * sscanf is undefined!
 * NOTE: sscanf can correctly handle any space in the target string as
 * well as in the format specifier, therefore, no need to process those
 * lines that will be processed by sscanf.
 * Footnote: In fprintf(), the rvalue type promotions are expected. %f and 
 * %g actually correspond to parameters of type double. Thus in fprintf()
 * there is no difference between %f and %lf, or between %g and %lg. However, 
 * in sscanf() what is passed is a pointer to the variable so no rvalue type 
 * promotions occur or are expected. Thus %f and %lf are quite different in
 * sscanf, but the same in fprintf. Consequently, we need to use %g for 
 * double in fprintf and %lg for double in sscanf. It doesn't matter which
 * you use for fprintf because the fprintf library function treats them as
 * synonymous, but it's crucial to get it right for sscanf. 
 */
static int ReadCaseSettingData(Space *space, Time *time, Model *model, Partition *part)
{
    FILE *filePointer = fopen("artracfd.case", "r");
    if (NULL == filePointer) {
        FatalError("failed to open case data file: artracfd.case...");
    }
    /*
     * Read file line by line to get case setting data
     */
    char currentLine[200] = {'\0'}; /* store the current read line */
    int entryCount = 0; /* entry count */
    /* set format specifier according to the type of Real */
    char formatI[5] = "%lg"; /* default is double type */
    char formatIII[15] = "%lg, %lg, %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
        strncpy(formatIII, "%g, %g, %g", sizeof formatIII); /* float type */
    }
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "space begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, 
                    &(space->xMin), &(space->yMin), &(space->zMin)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, 
                    &(space->xMax), &(space->yMax), &(space->zMax)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d, %d, %d", 
                    &(space->nx), &(space->ny), &(space->nz)); 
            continue;
        }
        if (0 == strncmp(currentLine, "time begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->restart)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(time->end)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->stepN)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(time->numCFL)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->outputN)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->dataStreamer)); 
            continue;
        }
        if (0 == strncmp(currentLine, "numerical begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(model->scheme)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(model->averager)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(model->splitter)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->delta)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(model->fsi)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(model->layers)); 
            continue;
        }
        if (0 == strncmp(currentLine, "fluid begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refPr)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refMu)); 
            continue;
        }
        if (0 == strncmp(currentLine, "reference begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refLength)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refDensity)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refVelocity)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refTemperature)); 
            continue;
        }
        if (0 == strncmp(currentLine, "initialization begin", sizeof currentLine)) {
            ++entryCount;
            /* Read initial values to inner part 0 */
            ReadConsecutiveRealData(&filePointer, part->valueBC[0], 5);
            continue;
        }
        if (0 == strncmp(currentLine, "west boundary begin", sizeof currentLine)) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 1);
            continue;
        }
        if (0 == strncmp(currentLine, "east boundary begin", sizeof currentLine)) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 2);
            continue;
        }
        if (0 == strncmp(currentLine, "south boundary begin", sizeof currentLine)) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 3);
            continue;
        }
        if (0 == strncmp(currentLine, "north boundary begin", sizeof currentLine)) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 4);
            continue;
        }
        if (0 == strncmp(currentLine, "front boundary begin", sizeof currentLine)) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 5);
            continue;
        }
        if (0 == strncmp(currentLine, "back boundary begin", sizeof currentLine)) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 6);
            continue;
        }
        if (0 == strncmp(currentLine, "plane initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->tallyIC] = 1; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->tallyIC] + 0, 
                    part->valueIC[part->tallyIC] + 1, part->valueIC[part->tallyIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->tallyIC] + 3, 
                    part->valueIC[part->tallyIC] + 4, part->valueIC[part->tallyIC] + 5); 
            ReadConsecutiveRealData(&filePointer, part->valueIC[part->tallyIC] + ENTRYIC - 5, 5);
            ++part->tallyIC; /* regional initializer count and pointer */
            continue;
        }
        if (0 == strncmp(currentLine, "sphere initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->tallyIC] = 2; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->tallyIC] + 0, 
                    part->valueIC[part->tallyIC] + 1, part->valueIC[part->tallyIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, part->valueIC[part->tallyIC] + 3); 
            ReadConsecutiveRealData(&filePointer, part->valueIC[part->tallyIC] + ENTRYIC - 5, 5);
            ++part->tallyIC; /* regional initializer count and pointer */
            continue;
        }
        if (0 == strncmp(currentLine, "box initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->tallyIC] = 3; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->tallyIC] + 0, 
                    part->valueIC[part->tallyIC] + 1, part->valueIC[part->tallyIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->tallyIC] + 3, 
                    part->valueIC[part->tallyIC] + 4, part->valueIC[part->tallyIC] + 5); 
            ReadConsecutiveRealData(&filePointer, part->valueIC[part->tallyIC] + ENTRYIC - 5, 5);
            ++part->tallyIC; /* regional initializer count and pointer */
            continue;
        }
        if (0 == strncmp(currentLine, "cylinder initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->tallyIC] = 4; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->tallyIC] + 0, 
                    part->valueIC[part->tallyIC] + 1, part->valueIC[part->tallyIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->tallyIC] + 3, 
                    part->valueIC[part->tallyIC] + 4, part->valueIC[part->tallyIC] + 5); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, part->valueIC[part->tallyIC] + 6); 
            ReadConsecutiveRealData(&filePointer, part->valueIC[part->tallyIC] + ENTRYIC - 5, 5);
            ++part->tallyIC; /* regional initializer count and pointer */
            continue;
        }
        if (0 == strncmp(currentLine, "probe control begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            fgets(currentLine, sizeof currentLine, filePointer);
            /* output time control of probes */
            sscanf(currentLine, "%d", &(time->outputProbe)); 
            continue;
        }
        if (0 == strncmp(currentLine, "probe begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, time->probe[time->tallyProbe] + 0, 
                    time->probe[time->tallyProbe] + 1, time->probe[time->tallyProbe] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, time->probe[time->tallyProbe] + 3, 
                    time->probe[time->tallyProbe] + 4, time->probe[time->tallyProbe] + 5); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, time->probe[time->tallyProbe] + 6);  /* resolution */
            ++time->tallyProbe; /* probe count and pointer */
            continue;
        }
    }
    fclose(filePointer); /* close current opened file */
    /*
     * Check missing information section in configuration
     */
    if (12 != entryCount) {
        FatalError("missing or repeated necessary information section");
    }
    return 0;
}
/*
 * Boundary type ID:
 * 0: interior region (default value)
 * 1: inflow
 * 2: outflow
 * 3: slip wall
 * 4: noslip wall
 * 5: periodic
 * -5: periodic pair
 */
static int ReadBoundaryData(FILE **filePointerPointer, Partition *part, const int partID)
{
    if (-5 == part->typeBC[partID]) { /* already set as periodic pair */
        return 0;
    }
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    char currentLine[200] = {'\0'}; /* store the current read line */
    char formatI[5] = "%lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    CommandLineProcessor(currentLine); /* process current line */
    if (0 == strncmp(currentLine, "inflow", sizeof currentLine)) {
        part->typeBC[partID] = 1;
        ReadConsecutiveRealData(&filePointer, part->valueBC[partID], 5);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "outflow", sizeof currentLine)) {
        part->typeBC[partID] = 2;
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "slip wall", sizeof currentLine)) {
        part->typeBC[partID] = 3;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[partID][5])); 
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "noslip wall", sizeof currentLine)) {
        part->typeBC[partID] = 4;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[partID][5])); 
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "periodic", sizeof currentLine)) {
        /* only need to set id and its periodic pair */
        part->typeBC[partID] = 5;
        if ((1 == partID) || (3 == partID) || (5 == partID)) {
            part->typeBC[partID+1] = -5;
        } else {
            part->typeBC[partID-1] = -5;
        }
        *filePointerPointer = filePointer;
        return 0;
    }
    FatalError("unidentified boundary type...");
    return 0;
}
/*
 * Read n consecutive real data entries into the memory pointed by address from
 * file pointed by the file pointer, and update the file pointer.
 */
static int ReadConsecutiveRealData(FILE **filePointerPointer, Real *address, const int entryN)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    char currentLine[200] = {'\0'}; /* store the current read line */
    char formatI[5] = "%lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
    }
    for (int n = 0; n < entryN; ++n) {
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, address + n); 
    }
    *filePointerPointer = filePointer; /* updated file pointer */
    return 0;
}
static int WriteBoundaryData(FILE **filePointerPointer, const Partition *part, const int partID)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    if (1 == part->typeBC[partID]) {
        fprintf(filePointer, "boundary type: inflow\n"); 
        fprintf(filePointer, "density: %.6g\n", part->valueBC[partID][0]);
        fprintf(filePointer, "x velocity: %.6g\n", part->valueBC[partID][1]);
        fprintf(filePointer, "y velocity: %.6g\n", part->valueBC[partID][2]);
        fprintf(filePointer, "z velocity: %.6g\n", part->valueBC[partID][3]);
        fprintf(filePointer, "pressure: %.6g\n", part->valueBC[partID][4]);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (2 == part->typeBC[partID]) {
        fprintf(filePointer, "boundary type: outflow\n"); 
        *filePointerPointer = filePointer;
        return 0;
    }
    if (3 == part->typeBC[partID]) {
        fprintf(filePointer, "boundary type: slip wall\n"); 
        fprintf(filePointer, "temperature: %.6g\n", part->valueBC[partID][5]);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (4 == part->typeBC[partID]) {
        fprintf(filePointer, "boundary type: noslip wall\n"); 
        fprintf(filePointer, "temperature: %.6g\n", part->valueBC[partID][5]);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (5 == part->typeBC[partID]) {
        fprintf(filePointer, "boundary type: periodic, primary pair with translation\n"); 
        *filePointerPointer = filePointer;
        return 0;
    }
    if (-5 == part->typeBC[partID]) {
        fprintf(filePointer, "boundary type: periodic, auxiliary pair with zero gradient flow\n"); 
        *filePointerPointer = filePointer;
        return 0;
    }
    FatalError("unidentified boundary type...");
    return 0;
}
static int WriteRegionalInitializerData(FILE **filePointerPointer, const Partition *part, const int n)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    if (1 == part->typeIC[n]) {
        fprintf(filePointer, "regional initialization: plane\n"); 
        fprintf(filePointer, "plane point x, y, z: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
        fprintf(filePointer, "plane normal nx, ny, nz: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][3], part->valueIC[n][4], part->valueIC[n][5]);
    }
    if (2 == part->typeIC[n]) {
        fprintf(filePointer, "regional initialization: sphere\n"); 
        fprintf(filePointer, "center point x, y, z: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
        fprintf(filePointer, "radius: %.6g\n", part->valueIC[n][3]);
    }
    if (3 == part->typeIC[n]) {
        fprintf(filePointer, "regional initialization: box\n"); 
        fprintf(filePointer, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
        fprintf(filePointer, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][3], part->valueIC[n][4], part->valueIC[n][5]);
    }
    if (4 == part->typeIC[n]) {
        fprintf(filePointer, "regional initialization: cylinder\n"); 
        fprintf(filePointer, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
        fprintf(filePointer, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][3], part->valueIC[n][4], part->valueIC[n][5]);
        fprintf(filePointer, "radius: %.6g\n", part->valueIC[n][6]);
    }
    fprintf(filePointer, "density: %.6g\n", part->valueIC[n][ENTRYIC-5]);
    fprintf(filePointer, "x velocity: %.6g\n", part->valueIC[n][ENTRYIC-4]);
    fprintf(filePointer, "y velocity: %.6g\n", part->valueIC[n][ENTRYIC-3]);
    fprintf(filePointer, "z velocity: %.6g\n", part->valueIC[n][ENTRYIC-2]);
    fprintf(filePointer, "pressure: %.6g\n", part->valueIC[n][ENTRYIC-1]);
    *filePointerPointer = filePointer;
    return 0;
}
/*
 * This function outputs the case setting data to a file for verification.
 */
static int WriteVerifyData(const Space *space, const Time *time, const Model *model, const Partition *part)
{
    ShowInformation("  Data outputted into artracfd.verify...");
    FILE *filePointer = fopen("artracfd.verify", "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to file: artracfd.verify");
    }
    /* output information to file */
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#                     Case Conformation for ArtraCFD                          -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Space Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "domain xmin, ymin, zmin: %.6g, %.6g, %.6g\n", space->xMin, space->yMin, space->zMin); 
    fprintf(filePointer, "domain xmax, ymax, zmax: %.6g, %.6g, %.6g\n", space->xMax, space->yMax, space->zMax); 
    fprintf(filePointer, "x, y, z mesh number: %d, %d, %d\n", space->nx, space->ny, space->nz); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Time Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "restart flag: %d\n", time->restart); 
    fprintf(filePointer, "total evolution time: %.6g\n", time->end); 
    fprintf(filePointer, "maximum number of steps: %d\n", time->stepN); 
    fprintf(filePointer, "CFL condition number: %.6g\n", time->numCFL); 
    fprintf(filePointer, "exporting data times: %d\n", time->outputN); 
    fprintf(filePointer, "data streamer id: %d\n", time->dataStreamer); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Numerical Method <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "spatial scheme: %d\n", model->scheme);
    fprintf(filePointer, "average method: %d\n", model->averager);
    fprintf(filePointer, "flux splitting method: %d\n", model->splitter);
    fprintf(filePointer, "Harten's numerical dissipation coefficient: %.6g\n", model->delta); 
    fprintf(filePointer, "fluid solid interaction: %d\n", model->fsi);
    fprintf(filePointer, "layers of ghost nodes using method of image: %d\n", model->layers);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Fluid and Flow Properties <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "Prandtl number: %.6g\n", model->refPr); 
    fprintf(filePointer, "modify coefficient of dynamic viscosity: %.6g\n", model->refMu); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Reference Values  <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "length: %.6g\n", model->refLength); 
    fprintf(filePointer, "density: %.6g\n", model->refDensity); 
    fprintf(filePointer, "velocity: %.6g\n", model->refVelocity); 
    fprintf(filePointer, "temperature: %.6g\n", model->refTemperature); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                            >> NOTE <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# Values in following parts are relative to reference values. Hence, they need\n");
    fprintf(filePointer, "# to be normalized by the given reference values. Like pressure should be\n");
    fprintf(filePointer, "# normalized by reference density times reference velocity square.\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                     >> Flow Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "density: %.6g\n", part->valueBC[0][0]);
    fprintf(filePointer, "x velocity: %.6g\n", part->valueBC[0][1]);
    fprintf(filePointer, "y velocity: %.6g\n", part->valueBC[0][2]);
    fprintf(filePointer, "z velocity: %.6g\n", part->valueBC[0][3]);
    fprintf(filePointer, "pressure: %.6g\n", part->valueBC[0][4]);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                     >> Boundary Condition <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian West\n"); 
    WriteBoundaryData(&filePointer, part, 1);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian East\n"); 
    WriteBoundaryData(&filePointer, part, 2);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian South\n"); 
    WriteBoundaryData(&filePointer, part, 3);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian North\n"); 
    WriteBoundaryData(&filePointer, part, 4);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Front\n"); 
    WriteBoundaryData(&filePointer, part, 5);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Back\n"); 
    WriteBoundaryData(&filePointer, part, 6);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                  >> Regional Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    for (int n = 0; n < part->tallyIC; ++n) {
        fprintf(filePointer, "#\n");
        WriteRegionalInitializerData(&filePointer, part, n);
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Field Data Probes <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    for (int n = 0; n < time->tallyProbe; ++n) {
        if (0 == n) {
            fprintf(filePointer, "total number of times of exporting probe data: %d\n", time->outputProbe);
        }
        fprintf(filePointer, "#\n");
        fprintf(filePointer, "probe %d\n", (n + 1));
        fprintf(filePointer, "x, y, z of the first end point: %.6g, %.6g, %.6g\n", 
                time->probe[n][0], time->probe[n][1], time->probe[n][2]);
        fprintf(filePointer, "x, y, z of the second end point: %.6g, %.6g, %.6g\n", 
                time->probe[n][3], time->probe[n][4], time->probe[n][5]);
        fprintf(filePointer, "number of points on line: %.6g\n", time->probe[n][6]);
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * This function do some parameter checking
 */
static int CheckCaseSettingData(const Space *space, const Time *time, const Model *model, const Partition *part)
{
    ShowInformation("  Preliminary case data checking ...");
    /* space */
    if ((0 >= (space->xMax - space->xMin)) || (0 >= (space->yMax - space->yMin)) ||
            (0 >= (space->zMax - space->zMin))) {
        FatalError("wrong domian region values in case settings");
    }
    if ((1 > space->nz) || (1 > space->ny) || (1 > space->nx)) {
        FatalError("too small mesh values in case settings");
    }
    /* time */
    if ((0 > time->restart) || (1 < time->restart)|| (0 >= time->end)
            || (0 >= time->numCFL ) || (1 > time->outputN)) {
        FatalError("wrong values in time section of case settings");
    }
    /* numerical method */
    if ((0 > model->scheme) || (1 < model->scheme) || 
            (0 > model->averager) || (1 < model->averager) ||
            (0 > model->splitter) || (1 < model->splitter) ||
            (0 > model->fsi) || (1 < model->fsi) ||
            (0 > model->delta)) {
        FatalError("wrong values in numerical method of case settings");
    }
    /* fluid and flow */
    if ((0 >= model->refPr) || (0 > model->refMu)) {
        FatalError("wrong values in fluid and flow section of case settings");
    }
    /* reference */
    if ((0 >= model->refLength) || (0 >= model->refDensity) || 
            (0 >= model->refVelocity) || (0 >= model->refTemperature)) {
        FatalError("wrong values in reference section of case settings");
    }
    /* initialization */
    if ((0 > part->valueBC[0][0]) || (0 > part->valueBC[0][4])) {
        FatalError("wrong values in initialization section of case settings");
    }
    return 0;
}
/* a good practice: end file with a newline */

