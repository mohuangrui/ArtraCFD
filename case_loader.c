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
static int ReadCaseSettingData(Time *, Space *, Model *);
static int ReadGeometrySettingData(Geometry *);
static int ReadBoundaryData(FILE **, Space *, const int);
static int ReadConsecutiveRealData(FILE **, Real *, const int);
static int WriteBoundaryData(FILE **, const Space *, const int);
static int WriteInitializerData(FILE **, const Space *, const int);
static int WriteVerifyData(const Time *, const Space *, const Model *);
static int CheckCaseSettingData(const Time *, const Space *, const Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int LoadCaseSettingData(Time *time, Space *space, Model *model)
{
    ReadCaseSettingData(time, space, model);
    ReadGeometrySettingData(&(space->geo));
    WriteVerifyData(time, space, model);
    CheckCaseSettingData(time, space, model);
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
static int ReadCaseSettingData(Time *time, Space *space, Model *model)
{
    Partition *part = &(space->part);
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
            sscanf(currentLine, formatIII, &(part->domain[X][MIN]),
                    &(part->domain[Y][MIN]), &(part->domain[Z][MIN])); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, &(part->domain[X][MAX]),
                    &(part->domain[Y][MAX]), &(part->domain[Z][MAX])); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d, %d, %d", 
                    &(part->m[X]), &(part->m[Y]), &(part->m[Z])); 
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
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(model->gState)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, &(model->g[X]),
                    &(model->g[Y]), &(model->g[Z])); 
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
            part->countIC = 0; /* enforce global initialization first */
            part->typeIC[part->countIC] = ICGLOBAL; /* IC type id */
            ReadConsecutiveRealData(&filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
            continue;
        }
        if (0 == strncmp(currentLine, "west boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, PWB);
            continue;
        }
        if (0 == strncmp(currentLine, "east boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, PEB);
            continue;
        }
        if (0 == strncmp(currentLine, "south boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, PSB);
            continue;
        }
        if (0 == strncmp(currentLine, "north boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, PNB);
            continue;
        }
        if (0 == strncmp(currentLine, "front boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, PFB);
            continue;
        }
        if (0 == strncmp(currentLine, "back boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(&filePointer, space, PBB);
            continue;
        }
        if (0 == strncmp(currentLine, "plane initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->countIC] = ICPLANE; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 0, 
                    part->valueIC[part->countIC] + 1, part->valueIC[part->countIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 3, 
                    part->valueIC[part->countIC] + 4, part->valueIC[part->countIC] + 5); 
            ReadConsecutiveRealData(&filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
            continue;
        }
        if (0 == strncmp(currentLine, "sphere initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->countIC] = ICSPHERE; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 0, 
                    part->valueIC[part->countIC] + 1, part->valueIC[part->countIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, part->valueIC[part->countIC] + 6); 
            ReadConsecutiveRealData(&filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
            continue;
        }
        if (0 == strncmp(currentLine, "box initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->countIC] = ICBOX; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 0, 
                    part->valueIC[part->countIC] + 1, part->valueIC[part->countIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 3, 
                    part->valueIC[part->countIC] + 4, part->valueIC[part->countIC] + 5); 
            ReadConsecutiveRealData(&filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
            continue;
        }
        if (0 == strncmp(currentLine, "cylinder initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->countIC] = ICCYLINDER; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 0, 
                    part->valueIC[part->countIC] + 1, part->valueIC[part->countIC] + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 3, 
                    part->valueIC[part->countIC] + 4, part->valueIC[part->countIC] + 5); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, part->valueIC[part->countIC] + 6); 
            ReadConsecutiveRealData(&filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
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
static int ReadGeometrySettingData(Geometry *geo)
{
    FILE *filePointer = fopen("artracfd.geo", "r");
    if (NULL == filePointer) {
        FatalError("failed to open file: artracfd.geo...");
    }
    /* read and process file line by line */
    String currentLine = {'\0'}; /* store the current read line */
    int entryCount = 0; /* entry count */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "count begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(geo->sphereN));
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(geo->stlN));
            if (0 >= geo->sphereN) {
                geo->sphereN = 0;
            }
            if (0 >= geo->stlN) {
                geo->stlN = 0;
            }
            geo->totalN = geo->sphereN + geo->stlN;
            break;
        }
    }
    fclose(filePointer); /* close current opened file */
    /* Check missing information section in configuration */
    if (1 != entryCount) {
        FatalError("missing necessary information section");
    }
    return 0;
}
static int ReadBoundaryData(FILE **filePointerPointer, Space *space, const int boundary)
{
    Partition *part = &(space->part);
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    String currentLine = {'\0'}; /* store the current read line */
    char formatI[5] = "%lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    CommandLineProcessor(currentLine); /* process current line */
    if (0 == strncmp(currentLine, "inflow", sizeof currentLine)) {
        part->typeBC[boundary] = INFLOW;
        ReadConsecutiveRealData(&filePointer, part->valueBC[boundary], VARBC);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "outflow", sizeof currentLine)) {
        part->typeBC[boundary] = OUTFLOW;
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "slip wall", sizeof currentLine)) {
        part->typeBC[boundary] = SLIPWALL;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[boundary][ENTRYBC-1]));
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "noslip wall", sizeof currentLine)) {
        part->typeBC[boundary] = NOSLIPWALL;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[boundary][ENTRYBC-1]));
        *filePointerPointer = filePointer;
        return 0;
    }
    if (0 == strncmp(currentLine, "periodic", sizeof currentLine)) {
        part->typeBC[boundary] = PERIODIC;
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
    const Partition *part = &(space->part);
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    if (INFLOW == part->typeBC[boundary]) {
        fprintf(filePointer, "boundary type: inflow\n"); 
        fprintf(filePointer, "density: %.6g\n", part->valueBC[boundary][0]);
        fprintf(filePointer, "x velocity: %.6g\n", part->valueBC[boundary][1]);
        fprintf(filePointer, "y velocity: %.6g\n", part->valueBC[boundary][2]);
        fprintf(filePointer, "z velocity: %.6g\n", part->valueBC[boundary][3]);
        fprintf(filePointer, "pressure: %.6g\n", part->valueBC[boundary][4]);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (OUTFLOW == part->typeBC[boundary]) {
        fprintf(filePointer, "boundary type: outflow\n"); 
        *filePointerPointer = filePointer;
        return 0;
    }
    if (SLIPWALL == part->typeBC[boundary]) {
        fprintf(filePointer, "boundary type: slip wall\n"); 
        fprintf(filePointer, "temperature: %.6g\n", part->valueBC[boundary][ENTRYBC-1]);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (NOSLIPWALL == part->typeBC[boundary]) {
        fprintf(filePointer, "boundary type: noslip wall\n"); 
        fprintf(filePointer, "temperature: %.6g\n", part->valueBC[boundary][ENTRYBC-1]);
        *filePointerPointer = filePointer;
        return 0;
    }
    if (PERIODIC == part->typeBC[boundary]) {
        fprintf(filePointer, "boundary type: periodic\n"); 
        *filePointerPointer = filePointer;
        return 0;
    }
    FatalError("unidentified boundary type...");
    return 0;
}
static int WriteInitializerData(FILE **filePointerPointer, const Space *space, const int n)
{
    const Partition *part = &(space->part);
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    if (ICGLOBAL == part->typeIC[n]) {
        ; /* no extra information */
    }
    if (ICPLANE == part->typeIC[n]) {
        fprintf(filePointer, "regional initialization: plane\n"); 
        fprintf(filePointer, "plane point x, y, z: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
        fprintf(filePointer, "plane normal nx, ny, nz: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][3], part->valueIC[n][4], part->valueIC[n][5]);
    }
    if (ICSPHERE == part->typeIC[n]) {
        fprintf(filePointer, "regional initialization: sphere\n"); 
        fprintf(filePointer, "center point x, y, z: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
        fprintf(filePointer, "radius: %.6g\n", part->valueIC[n][6]);
    }
    if (ICBOX == part->typeIC[n]) {
        fprintf(filePointer, "regional initialization: box\n"); 
        fprintf(filePointer, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
        fprintf(filePointer, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n", 
                part->valueIC[n][3], part->valueIC[n][4], part->valueIC[n][5]);
    }
    if (ICCYLINDER == part->typeIC[n]) {
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
static int WriteVerifyData(const Time *time, const Space *space, const Model *model)
{
    const Partition *part = &(space->part);
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
    fprintf(filePointer, "domain xmin, ymin, zmin: %.6g, %.6g, %.6g\n", 
            part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]); 
    fprintf(filePointer, "domain xmax, ymax, zmax: %.6g, %.6g, %.6g\n", 
            part->domain[X][MAX], part->domain[Y][MAX], part->domain[Z][MAX]); 
    fprintf(filePointer, "x, y, z mesh number: %d, %d, %d\n", part->m[X], part->m[Y], part->m[Z]); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Time Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "restart number tag: %d\n", time->restart); 
    fprintf(filePointer, "total evolution time: %.6g\n", time->end); 
    fprintf(filePointer, "CFL condition number: %.6g\n", time->numCFL); 
    fprintf(filePointer, "maximum number of steps: %d\n", time->stepN); 
    fprintf(filePointer, "exporting data times: %d\n", time->outputN); 
    fprintf(filePointer, "data streamer: %d\n", time->dataStreamer); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Numerical Method <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "spatial scheme: %d\n", model->scheme);
    fprintf(filePointer, "average method: %d\n", model->averager);
    fprintf(filePointer, "flux splitting method: %d\n", model->splitter);
    fprintf(filePointer, "phase interaction: %d\n", model->fsi);
    fprintf(filePointer, "interfacial layers using reconstruction: %d\n", model->layers);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                       >> Material Properties <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "material: %d\n", model->matID); 
    fprintf(filePointer, "viscous level: %.6g\n", model->refMu); 
    fprintf(filePointer, "gravity state: %d\n", model->gState); 
    fprintf(filePointer, "gravity vector: %.6g, %.6g, %.6g\n", model->g[X], model->g[Y], model->g[Z]); 
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
    WriteBoundaryData(&filePointer, space, PWB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian East\n"); 
    WriteBoundaryData(&filePointer, space, PEB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian South\n"); 
    WriteBoundaryData(&filePointer, space, PSB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian North\n"); 
    WriteBoundaryData(&filePointer, space, PNB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Front\n"); 
    WriteBoundaryData(&filePointer, space, PFB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Back\n"); 
    WriteBoundaryData(&filePointer, space, PBB);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                  >> Regional Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    for (int n = 1; n < part->countIC; ++n) {
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
static int CheckCaseSettingData(const Time *time, const Space *space, const Model *model)
{
    const Partition *part = &(space->part);
    /* space */
    if ((0.0 >= (part->domain[X][MAX] - part->domain[X][MIN])) ||
            (0.0 >= (part->domain[Y][MAX] - part->domain[Y][MIN])) ||
            (0.0 >= (part->domain[Z][MAX] - part->domain[Z][MIN]))) {
        FatalError("wrong domian region values in case settings");
    }
    if ((1 > part->m[X]) || (1 > part->m[Y]) || (1 > part->m[Z])) {
        FatalError("too small mesh values in case settings");
    }
    /* time */
    if ((0 > time->restart) || (0.0 >= time->end)
            || (0.0 >= time->numCFL) || (1 > time->outputN)) {
        FatalError("wrong values in time section of case settings");
    }
    /* numerical method */
    if ((0 > model->scheme) || (0 > model->averager) ||
            (0 > model->splitter) || (0 > model->fsi)) {
        FatalError("wrong values in numerical method of case settings");
    }
    /* material */
    if ((0 > model->matID)) {
        FatalError("wrong values in material section of case settings");
    }
    /* reference */
    if ((0.0 >= model->refL) || (0.0 >= model->refRho) || 
            (0.0 >= model->refV) || (0.0 >= model->refT)) {
        FatalError("wrong values in reference section of case settings");
    }
    return 0;
}
/* a good practice: end file with a newline */

