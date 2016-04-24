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
#include "case_loader.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ReadCaseSettingData(Space *, Time *, Model *);
static int ReadBoundaryData(FILE **, Space *, const int);
static int ReadConsecutiveRealData(FILE **, Real *, const int);
static int WriteBoundaryData(FILE **, const Space *, const int);
static int WriteInitializerData(FILE **, const Space *, const int);
static int WriteVerifyData(const Space *, const Time *, const Model *);
static int CheckCaseSettingData(const Space *, const Time *, const Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int LoadCaseSettingData(Space *space, Time *time, Model *model)
{
    ShowInformation("Loading case setting data ...");
    ReadCaseSettingData(space, time, model);
    WriteVerifyData(space, time, model);
    CheckCaseSettingData(space, time, model);
    ShowInformation("Session End");
    return 0;
}
/*
 * This function read the case settings from the case file.
 * The key is to read and process file line by line. Use "begin end" 
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
 * float and %lg for double in sscanf. It doesn't matter which
 * you use for fprintf because the fprintf library function treats them as
 * synonymous, but it's crucial to get it right for sscanf. 
 */
static int ReadCaseSettingData(Space *space, Time *time, Model *model)
{
    FILE *filePointer = fopen("artracfd.case", "r");
    if (NULL == filePointer) {
        FatalError("failed to open case data file: artracfd.case...");
    }
    /*
     * Read file line by line to get case setting data
     */
    String currentLine = {'\0'}; /* store the current read line */
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
                    &(space->part.domain[X][MIN]), &(space->part.domain[Y][MIN]), &(space->part.domain[Z][MIN])); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, 
                    &(space->part.domain[X][MAX]), &(space->part.domain[Y][MAX]), &(space->part.domain[Z][MAX])); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d, %d, %d", 
                    &(space->part.m[X]), &(space->part.m[Y]), &(space->part.m[Z])); 
            continue;
        }
        if (0 == strncmp(currentLine, "time begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->restart)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(time->end)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(time->numCFL)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->stepN)); 
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
            sscanf(currentLine, "%d", &(model->fsi)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(model->layers)); 
            continue;
        }
        if (0 == strncmp(currentLine, "material begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(model->matID)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refMu)); 
            continue;
        }
        if (0 == strncmp(currentLine, "reference begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refL)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refRho)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refV)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(model->refT)); 
            continue;
        }
        if (0 == strncmp(currentLine, "initialization begin", sizeof currentLine)) {
            ++entryCount;
            space->part.typeIC[space->part.countIC] = ICGLOBAL; /* IC type id */
            ReadConsecutiveRealData(&filePointer, space->part.valueIC[space->part.countIC] + ENTRYIC - VARIC, VARIC);
            ++space->part.countIC; /* initializer count and pointer */
            continue;
        }
        if (0 == strncmp(currentLine, "west boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, BCWEST);
            continue;
        }
        if (0 == strncmp(currentLine, "east boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, BCEAST);
            continue;
        }
        if (0 == strncmp(currentLine, "south boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, BCSOUTH);
            continue;
        }
        if (0 == strncmp(currentLine, "north boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, BCNORTH);
            continue;
        }
        if (0 == strncmp(currentLine, "front boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, BCFRONT);
            continue;
        }
        if (0 == strncmp(currentLine, "back boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, BCBACK);
            continue;
        }
        if (0 == strncmp(currentLine, "plane initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            space->part.typeIC[space->part.countIC] = ICPLANE; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, space->part.valueIC[space->part.countIC] + 0, 
                    space->part.valueIC[space->part.countIC] + 1, space->part.valueIC[space->part.countIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, space->part.valueIC[space->part.countIC] + 3, 
                    space->part.valueIC[space->part.countIC] + 4, space->part.valueIC[space->part.countIC] + 5); 
            ReadConsecutiveRealData(&filePointer, space->part.valueIC[space->part.countIC] + ENTRYIC - VARIC, VARIC);
            ++space->part.countIC; /* initializer count and pointer */
            continue;
        }
        if (0 == strncmp(currentLine, "sphere initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            space->part.typeIC[space->part.countIC] = ICSPHERE; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, space->part.valueIC[space->part.countIC] + 0, 
                    space->part.valueIC[space->part.countIC] + 1, space->part.valueIC[space->part.countIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, space->part.valueIC[space->part.countIC] + 3); 
            ReadConsecutiveRealData(&filePointer, space->part.valueIC[space->part.countIC] + ENTRYIC - VARIC, VARIC);
            ++space->part.countIC; /* initializer count and pointer */
            continue;
        }
        if (0 == strncmp(currentLine, "box initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            space->part.typeIC[space->part.countIC] = ICBOX; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, space->part.valueIC[space->part.countIC] + 0, 
                    space->part.valueIC[space->part.countIC] + 1, space->part.valueIC[space->part.countIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, space->part.valueIC[space->part.countIC] + 3, 
                    space->part.valueIC[space->part.countIC] + 4, space->part.valueIC[space->part.countIC] + 5); 
            ReadConsecutiveRealData(&filePointer, space->part.valueIC[space->part.countIC] + ENTRYIC - VARIC, VARIC);
            ++space->part.countIC; /* initializer count and pointer */
            continue;
        }
        if (0 == strncmp(currentLine, "cylinder initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            space->part.typeIC[space->part.countIC] = ICCYLINDER; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, space->part.valueIC[space->part.countIC] + 0, 
                    space->part.valueIC[space->part.countIC] + 1, space->part.valueIC[space->part.countIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, space->part.valueIC[space->part.countIC] + 3, 
                    space->part.valueIC[space->part.countIC] + 4, space->part.valueIC[space->part.countIC] + 5); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, space->part.valueIC[space->part.countIC] + 6); 
            ReadConsecutiveRealData(&filePointer, space->part.valueIC[space->part.countIC] + ENTRYIC - VARIC, VARIC);
            ++space->part.countIC; /* initializer count and pointer */
            continue;
        }
        if (0 == strncmp(currentLine, "probe control begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            fgets(currentLine, sizeof currentLine, filePointer);
            /* output time control of probes */
            sscanf(currentLine, "%d", &(time->outputNProbe)); 
            continue;
        }
        if (0 == strncmp(currentLine, "probe begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, time->probe[time->probeN] + 0, 
                    time->probe[time->probeN] + 1, time->probe[time->probeN] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, time->probe[time->probeN] + 3, 
                    time->probe[time->probeN] + 4, time->probe[time->probeN] + 5); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, time->probe[time->probeN] + 6);  /* resolution */
            ++time->probeN; /* probe count and pointer */
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
static int ReadBoundaryData(FILE **filePointerPointer, Space *space, const int boundary)
{
    if (PERIODICPAIR == space->part.typeBC[boundary]) { /* already set as periodic pair */
        return 0;
    }
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    String currentLine = {'\0'}; /* store the current read line */
    char formatI[5] = "%lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    CommandLineProcessor(currentLine); /* process current line */
    if (0 == strncmp(currentLine, "inflow", sizeof currentLine)) {
        space->part.typeBC[boundary] = INFLOW;
        ReadConsecutiveRealData(&filePointer, space->part.valueBC[boundary], VARBC);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "outflow", sizeof currentLine)) {
        space->part.typeBC[boundary] = OUTFLOW;
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "slip wall", sizeof currentLine)) {
        space->part.typeBC[boundary] = SLIPWALL;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(space->part.valueBC[boundary][ENTRYBC-1]));
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "noslip wall", sizeof currentLine)) {
        space->part.typeBC[boundary] = NOSLIPWALL;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(space->part.valueBC[boundary][ENTRYBC-1]));
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "periodic", sizeof currentLine)) {
        /* only need to set id and its periodic pair */
        space->part.typeBC[boundary] = PERIODIC;
        if ((BCWEST == boundary) || (BCSOUTH == boundary) || (BCFRONT == boundary)) {
            space->part.typeBC[boundary+1] = PERIODICPAIR;
        } else {
            space->part.typeBC[boundary-1] = PERIODICPAIR;
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
    String currentLine = {'\0'}; /* store the current read line */
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
static int WriteBoundaryData(FILE **filePointerPointer, const Space *space, const int boundary)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    if (INFLOW == space->part.typeBC[boundary]) {
        fprintf(filePointer, "boundary type: inflow\n"); 
        fprintf(filePointer, "density: %.6g\n", space->part.valueBC[boundary][0]);
        fprintf(filePointer, "x velocity: %.6g\n", space->part.valueBC[boundary][1]);
        fprintf(filePointer, "y velocity: %.6g\n", space->part.valueBC[boundary][2]);
        fprintf(filePointer, "z velocity: %.6g\n", space->part.valueBC[boundary][3]);
        fprintf(filePointer, "pressure: %.6g\n", space->part.valueBC[boundary][4]);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (OUTFLOW == space->part.typeBC[boundary]) {
        fprintf(filePointer, "boundary type: outflow\n"); 
        *filePointerPointer = filePointer;
        return 0;
    }
    if (SLIPWALL == space->part.typeBC[boundary]) {
        fprintf(filePointer, "boundary type: slip wall\n"); 
        fprintf(filePointer, "temperature: %.6g\n", space->part.valueBC[boundary][ENTRYBC-1]);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (NOSLIPWALL == space->part.typeBC[boundary]) {
        fprintf(filePointer, "boundary type: noslip wall\n"); 
        fprintf(filePointer, "temperature: %.6g\n", space->part.valueBC[boundary][ENTRYBC-1]);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (PERIODIC == space->part.typeBC[boundary]) {
        fprintf(filePointer, "boundary type: periodic, primary pair with translation\n"); 
        *filePointerPointer = filePointer;
        return 0;
    }
    if (PERIODICPAIR == space->part.typeBC[boundary]) {
        fprintf(filePointer, "boundary type: periodic, auxiliary pair\n"); 
        *filePointerPointer = filePointer;
        return 0;
    }
    FatalError("unidentified boundary type...");
    return 0;
}
static int WriteInitializerData(FILE **filePointerPointer, const Space *space, const int n)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    if (ICGLOBAL == space->part.typeIC[n]) {
        ; /* no extra information */
    }
    if (ICPLANE == space->part.typeIC[n]) {
        fprintf(filePointer, "regional initialization: plane\n"); 
        fprintf(filePointer, "plane point x, y, z: %.6g, %.6g, %.6g\n", 
                space->part.valueIC[n][0], space->part.valueIC[n][1], space->part.valueIC[n][2]);
        fprintf(filePointer, "plane normal nx, ny, nz: %.6g, %.6g, %.6g\n", 
                space->part.valueIC[n][3], space->part.valueIC[n][4], space->part.valueIC[n][5]);
    }
    if (ICSPHERE == space->part.typeIC[n]) {
        fprintf(filePointer, "regional initialization: sphere\n"); 
        fprintf(filePointer, "center point x, y, z: %.6g, %.6g, %.6g\n", 
                space->part.valueIC[n][0], space->part.valueIC[n][1], space->part.valueIC[n][2]);
        fprintf(filePointer, "radius: %.6g\n", space->part.valueIC[n][3]);
    }
    if (ICBOX == space->part.typeIC[n]) {
        fprintf(filePointer, "regional initialization: box\n"); 
        fprintf(filePointer, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n", 
                space->part.valueIC[n][0], space->part.valueIC[n][1], space->part.valueIC[n][2]);
        fprintf(filePointer, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n", 
                space->part.valueIC[n][3], space->part.valueIC[n][4], space->part.valueIC[n][5]);
    }
    if (ICCYLINDER == space->part.typeIC[n]) {
        fprintf(filePointer, "regional initialization: cylinder\n"); 
        fprintf(filePointer, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n", 
                space->part.valueIC[n][0], space->part.valueIC[n][1], space->part.valueIC[n][2]);
        fprintf(filePointer, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n", 
                space->part.valueIC[n][3], space->part.valueIC[n][4], space->part.valueIC[n][5]);
        fprintf(filePointer, "radius: %.6g\n", space->part.valueIC[n][6]);
    }
    fprintf(filePointer, "density: %.6g\n", space->part.valueIC[n][ENTRYIC-5]);
    fprintf(filePointer, "x velocity: %.6g\n", space->part.valueIC[n][ENTRYIC-4]);
    fprintf(filePointer, "y velocity: %.6g\n", space->part.valueIC[n][ENTRYIC-3]);
    fprintf(filePointer, "z velocity: %.6g\n", space->part.valueIC[n][ENTRYIC-2]);
    fprintf(filePointer, "pressure: %.6g\n", space->part.valueIC[n][ENTRYIC-1]);
    *filePointerPointer = filePointer;
    return 0;
}
/*
 * This function outputs the case setting data to a file for verification.
 */
static int WriteVerifyData(const Space *space, const Time *time, const Model *model)
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
    fprintf(filePointer, "domain xmin, ymin, zmin: %.6g, %.6g, %.6g\n", space->part.domain[X][MIN], space->part.domain[Y][MIN], space->part.domain[Z][MIN]); 
    fprintf(filePointer, "domain xmax, ymax, zmax: %.6g, %.6g, %.6g\n", space->part.domain[X][MAX], space->part.domain[Y][MAX], space->part.domain[Z][MAX]); 
    fprintf(filePointer, "x, y, z mesh number: %d, %d, %d\n", space->part.m[X], space->part.m[Y], space->part.m[Z]); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Time Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "restart number tag: %d\n", time->restart); 
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
    fprintf(filePointer, "material interaction: %d\n", model->fsi);
    fprintf(filePointer, "layers of ghost nodes using method of image: %d\n", model->layers);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                       >> Material Properties <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "material: %d\n", model->matID); 
    fprintf(filePointer, "viscous level: %.6g\n", model->refMu); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Reference Values  <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "length: %.6g\n", model->refL); 
    fprintf(filePointer, "density: %.6g\n", model->refRho); 
    fprintf(filePointer, "velocity: %.6g\n", model->refV); 
    fprintf(filePointer, "temperature: %.6g\n", model->refT); 
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
    WriteInitializerData(&filePointer, space, ICGLOBAL);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                     >> Boundary Condition <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian West\n"); 
    WriteBoundaryData(&filePointer, space, BCWEST);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian East\n"); 
    WriteBoundaryData(&filePointer, space, BCEAST);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian South\n"); 
    WriteBoundaryData(&filePointer, space, BCSOUTH);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian North\n"); 
    WriteBoundaryData(&filePointer, space, BCNORTH);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Front\n"); 
    WriteBoundaryData(&filePointer, space, BCFRONT);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Back\n"); 
    WriteBoundaryData(&filePointer, space, BCBACK);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                  >> Regional Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    for (int n = 1; n < space->part.countIC; ++n) {
        fprintf(filePointer, "#\n");
        WriteInitializerData(&filePointer, space, n);
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Field Data Probes <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    for (int n = 0; n < time->probeN; ++n) {
        if (0 == n) {
            fprintf(filePointer, "total number of times of exporting probe data: %d\n", time->outputNProbe);
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
static int CheckCaseSettingData(const Space *space, const Time *time, const Model *model)
{
    ShowInformation("  Preliminary case data checking ...");
    /* space */
    if ((0 >= (space->part.domain[X][MAX] - space->part.domain[X][MIN])) ||
            (0 >= (space->part.domain[Y][MAX] - space->part.domain[Y][MIN])) ||
            (0 >= (space->part.domain[Z][MAX] - space->part.domain[Z][MIN]))) {
        FatalError("wrong domian region values in case settings");
    }
    if ((1 > space->part.m[X]) || (1 > space->part.m[Y]) || (1 > space->part.m[Z])) {
        FatalError("too small mesh values in case settings");
    }
    /* time */
    if ((0 > time->restart) || (0 >= time->end)
            || (0 >= time->numCFL) || (1 > time->outputN)) {
        FatalError("wrong values in time section of case settings");
    }
    /* numerical method */
    if ((0 > model->scheme) || (1 < model->scheme) || 
            (0 > model->averager) || (1 < model->averager) ||
            (0 > model->splitter) || (1 < model->splitter) ||
            (0 > model->fsi) || (1 < model->fsi)) {
        FatalError("wrong values in numerical method of case settings");
    }
    /* material */
    if ((0 > model->matID) || (2 < model->matID) || (0 > model->refMu)) {
        FatalError("wrong values in material section of case settings");
    }
    /* reference */
    if ((0 >= model->refL) || (0 >= model->refRho) || 
            (0 >= model->refV) || (0 >= model->refT)) {
        FatalError("wrong values in reference section of case settings");
    }
    return 0;
}
/* a good practice: end file with a newline */

