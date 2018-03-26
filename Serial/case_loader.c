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
static int ReadBoundaryData(FILE *, Space *, const int);
static int ReadConsecutiveRealData(FILE *, Real *, const int);
static int WriteBoundaryData(FILE *, const Space *, const int);
static int WriteInitializerData(FILE *, const Space *, const int);
static int WriteVerifyData(const Time *, const Space *, const Model *);
static int CheckCaseSettingData(const Time *, const Space *, const Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int LoadCaseData(Time *time, Space *space, Model *model)
{
    ReadCaseSettingData(time, space, model);
    ReadGeometrySettingData(&(space->geo));
    WriteVerifyData(time, space, model);
    CheckCaseSettingData(time, space, model);
    return 0;
}
/*
 * This function reads the case settings from the case file.
 * Use "begin end" environment in the case file to control the reading. 
 * The function scanf is notorious for its poor end-of-line handling.
 * Instead, use fgets to read a line of input and sscanf to process it.
 * Note: use a large enough number when using fgets to ensure reading
 * a whole line at a time. fgets will get the entire line including
 * the newline character (\n).
 * Note: sscanf can correctly handle any space in the target string as
 * well as in the format specifier, therefore, no need to process those
 * lines that will be processed by sscanf.
 * Note: In fprintf(), the rvalue type promotions are expected. %f and 
 * %g actually correspond to parameters of type double. Thus in fprintf()
 * there is no difference between %f and %lf, or between %g and %lg. However, 
 * in sscanf() what is passed is a pointer to the variable so no rvalue type 
 * promotions occur or are expected. Thus %f and %lf are quite different in
 * sscanf, but the same in fprintf. Consequently, we need to use %g for 
 * float and %lg for double in sscanf. It doesn't matter which you use for 
 * fprintf because the fprintf library function treats them as synonymous, 
 * but it's crucial to get it right for sscanf. 
 */
static int ReadCaseSettingData(Time *time, Space *space, Model *model)
{
    Partition *part = &(space->part);
    FILE *filePointer = fopen("artracfd.case", "r");
    if (NULL == filePointer) {
        FatalError("failed to open file: artracfd.case...");
    }
    /*
     * Read file line by line to get case setting data
     */
    String currentLine = {'\0'}; /* store the current read line */
    int nscan = 0; /* read conversion count */
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
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, &(part->domain[X][MIN]),
                    &(part->domain[Y][MIN]), &(part->domain[Z][MIN])); 
            VerifyReadConversion(nscan, 3);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, &(part->domain[X][MAX]),
                    &(part->domain[Y][MAX]), &(part->domain[Z][MAX])); 
            VerifyReadConversion(nscan, 3);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d, %d, %d", 
                    &(part->m[X]), &(part->m[Y]), &(part->m[Z])); 
            VerifyReadConversion(nscan, 3);
            continue;
        }
        if (0 == strncmp(currentLine, "time begin", sizeof currentLine)) {
            ++entryCount;
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->restart)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatI, &(time->end)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatI, &(time->numCFL)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->stepN)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->writeN)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->dataStreamer)); 
            VerifyReadConversion(nscan, 1);
            continue;
        }
        if (0 == strncmp(currentLine, "numerical begin", sizeof currentLine)) {
            ++entryCount;
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(model->tScheme)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(model->sScheme)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(model->multidim)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(model->jacobMean)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(model->fluxSplit)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(model->fsi)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(model->ibmLayer)); 
            VerifyReadConversion(nscan, 1);
            continue;
        }
        if (0 == strncmp(currentLine, "material begin", sizeof currentLine)) {
            ++entryCount;
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(model->mid)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatI, &(model->refMu)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(model->gState)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, &(model->g[X]),
                    &(model->g[Y]), &(model->g[Z])); 
            VerifyReadConversion(nscan, 3);
            continue;
        }
        if (0 == strncmp(currentLine, "reference begin", sizeof currentLine)) {
            ++entryCount;
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatI, &(model->refL)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatI, &(model->refRho)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatI, &(model->refV)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatI, &(model->refT)); 
            VerifyReadConversion(nscan, 1);
            continue;
        }
        if (0 == strncmp(currentLine, "initialization begin", sizeof currentLine)) {
            ++entryCount;
            part->countIC = 0; /* enforce global initialization first */
            part->typeIC[part->countIC] = ICGLOBAL; /* IC type id */
            ReadConsecutiveRealData(filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
            continue;
        }
        if (0 == strncmp(currentLine, "west boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(filePointer, space, PWB);
            continue;
        }
        if (0 == strncmp(currentLine, "east boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(filePointer, space, PEB);
            continue;
        }
        if (0 == strncmp(currentLine, "south boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(filePointer, space, PSB);
            continue;
        }
        if (0 == strncmp(currentLine, "north boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(filePointer, space, PNB);
            continue;
        }
        if (0 == strncmp(currentLine, "front boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(filePointer, space, PFB);
            continue;
        }
        if (0 == strncmp(currentLine, "back boundary begin", sizeof currentLine)) {
            ++entryCount;
            ReadBoundaryData(filePointer, space, PBB);
            continue;
        }
        if (0 == strncmp(currentLine, "plane initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->countIC] = ICPLANE; /* IC type id */
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 0, 
                    part->valueIC[part->countIC] + 1, part->valueIC[part->countIC] + 2); 
            VerifyReadConversion(nscan, 3);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 3, 
                    part->valueIC[part->countIC] + 4, part->valueIC[part->countIC] + 5); 
            VerifyReadConversion(nscan, 3);
            ReadConsecutiveRealData(filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
            continue;
        }
        if (0 == strncmp(currentLine, "sphere initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->countIC] = ICSPHERE; /* IC type id */
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 0, 
                    part->valueIC[part->countIC] + 1, part->valueIC[part->countIC] + 2); 
            VerifyReadConversion(nscan, 3);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatI, part->valueIC[part->countIC] + 6); 
            VerifyReadConversion(nscan, 1);
            ReadConsecutiveRealData(filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
            continue;
        }
        if (0 == strncmp(currentLine, "box initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->countIC] = ICBOX; /* IC type id */
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 0, 
                    part->valueIC[part->countIC] + 1, part->valueIC[part->countIC] + 2); 
            VerifyReadConversion(nscan, 3);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 3, 
                    part->valueIC[part->countIC] + 4, part->valueIC[part->countIC] + 5); 
            VerifyReadConversion(nscan, 3);
            ReadConsecutiveRealData(filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
            continue;
        }
        if (0 == strncmp(currentLine, "cylinder initialization begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            part->typeIC[part->countIC] = ICCYLINDER; /* IC type id */
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 0, 
                    part->valueIC[part->countIC] + 1, part->valueIC[part->countIC] + 2); 
            VerifyReadConversion(nscan, 3);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatIII, part->valueIC[part->countIC] + 3, 
                    part->valueIC[part->countIC] + 4, part->valueIC[part->countIC] + 5); 
            VerifyReadConversion(nscan, 3);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, formatI, part->valueIC[part->countIC] + 6); 
            VerifyReadConversion(nscan, 1);
            ReadConsecutiveRealData(filePointer, part->valueIC[part->countIC] + ENTRYIC - VARIC, VARIC);
            ++part->countIC;
            continue;
        }
        if (0 == strncmp(currentLine, "probe count begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->pointProbeN)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->lineProbeN)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->curveProbeN)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->forceProbeN)); 
            VerifyReadConversion(nscan, 1);
            if (0 < time->pointProbeN) {
                time->pp = AssignStorage(time->pointProbeN * sizeof(*time->pp));
            }
            if (0 < time->lineProbeN) {
                time->lp = AssignStorage(time->lineProbeN * sizeof(*time->lp));
            }
            continue;
        }
        if (0 == strncmp(currentLine, "probe control begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->pointWriteN)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->lineWriteN)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->curveWriteN)); 
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(time->forceWriteN)); 
            VerifyReadConversion(nscan, 1);
            continue;
        }
        if (0 == strncmp(currentLine, "point probe begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            for (int n = 0; n < time->pointProbeN; ++n) {
                Fgets(currentLine, sizeof currentLine, filePointer);
                nscan = sscanf(currentLine, formatIII, time->pp[n] + 0, 
                        time->pp[n] + 1, time->pp[n] + 2);
                VerifyReadConversion(nscan, 3);
            }
            continue;
        }
        if (0 == strncmp(currentLine, "line probe begin", sizeof currentLine)) {
            /* optional entry do not increase entry count */
            for (int n = 0; n < time->lineProbeN; ++n) {
                Fgets(currentLine, sizeof currentLine, filePointer);
                nscan = sscanf(currentLine, formatIII, time->lp[n] + 0, 
                        time->lp[n] + 1, time->lp[n] + 2);
                VerifyReadConversion(nscan, 3);
                Fgets(currentLine, sizeof currentLine, filePointer);
                nscan = sscanf(currentLine, formatIII, time->lp[n] + 3, 
                        time->lp[n] + 4, time->lp[n] + 5);
                VerifyReadConversion(nscan, 3);
                Fgets(currentLine, sizeof currentLine, filePointer);
                nscan = sscanf(currentLine, formatI, time->lp[n] + 6);
                VerifyReadConversion(nscan, 1);
            }
            continue;
        }
    }
    fclose(filePointer); /* close current opened file */
    /*
     * Check missing information section in configuration
     */
    if (12 != entryCount) {
        FatalError("missing or repeated sections in case file");
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
    int nscan = 0; /* read conversion count */
    int entryCount = 0; /* entry count */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "count begin", sizeof currentLine)) {
            ++entryCount;
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(geo->sphN));
            VerifyReadConversion(nscan, 1);
            Fgets(currentLine, sizeof currentLine, filePointer);
            nscan = sscanf(currentLine, "%d", &(geo->stlN));
            VerifyReadConversion(nscan, 1);
            break;
        }
    }
    fclose(filePointer); /* close current opened file */
    /* Check missing information section in configuration */
    if (1 != entryCount) {
        FatalError("missing or repeated sections in case file");
    }
    return 0;
}
static int ReadBoundaryData(FILE *filePointer, Space *space, const int n)
{
    Partition *part = &(space->part);
    String currentLine = {'\0'}; /* store the current read line */
    int nscan = 0; /* read conversion count */
    char formatI[5] = "%lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
    }
    Fgets(currentLine, sizeof currentLine, filePointer);
    CommandLineProcessor(currentLine); /* process current line */
    if (0 == strncmp(currentLine, "inflow", sizeof currentLine)) {
        part->typeBC[n] = INFLOW;
        ReadConsecutiveRealData(filePointer, part->valueBC[n], VARBC);
        return 0;
    }
    if (0 == strncmp(currentLine, "outflow", sizeof currentLine)) {
        part->typeBC[n] = OUTFLOW;
        return 0;
    }
    if (0 == strncmp(currentLine, "slip wall", sizeof currentLine)) {
        part->typeBC[n] = SLIPWALL;
        Fgets(currentLine, sizeof currentLine, filePointer);
        nscan = sscanf(currentLine, formatI, &(part->valueBC[n][ENTRYBC-1]));
        VerifyReadConversion(nscan, 1);
        return 0;
    }
    if (0 == strncmp(currentLine, "noslip wall", sizeof currentLine)) {
        part->typeBC[n] = NOSLIPWALL;
        Fgets(currentLine, sizeof currentLine, filePointer);
        nscan = sscanf(currentLine, formatI, &(part->valueBC[n][ENTRYBC-1]));
        VerifyReadConversion(nscan, 1);
        return 0;
    }
    if (0 == strncmp(currentLine, "periodic", sizeof currentLine)) {
        part->typeBC[n] = PERIODIC;
        return 0;
    }
    FatalError("unidentified boundary type...");
    return 0;
}
/*
 * Read n consecutive real data entries into the memory pointed by address from
 * file pointed by the file pointer, and update the file pointer.
 */
static int ReadConsecutiveRealData(FILE *filePointer, Real *address, const int entryN)
{
    String currentLine = {'\0'}; /* store the current read line */
    int nscan = 0; /* read conversion count */
    char formatI[5] = "%lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
    }
    for (int n = 0; n < entryN; ++n) {
        Fgets(currentLine, sizeof currentLine, filePointer);
        nscan = sscanf(currentLine, formatI, address + n); 
        VerifyReadConversion(nscan, 1);
    }
    return 0;
}
static int WriteBoundaryData(FILE *filePointer, const Space *space, const int n)
{
    const Partition *part = &(space->part);
    switch (part->typeBC[n]) {
        case INFLOW:
            fprintf(filePointer, "boundary type: inflow\n"); 
            fprintf(filePointer, "density: %.6g\n", part->valueBC[n][0]);
            fprintf(filePointer, "x velocity: %.6g\n", part->valueBC[n][1]);
            fprintf(filePointer, "y velocity: %.6g\n", part->valueBC[n][2]);
            fprintf(filePointer, "z velocity: %.6g\n", part->valueBC[n][3]);
            fprintf(filePointer, "pressure: %.6g\n", part->valueBC[n][4]);
            break;
        case OUTFLOW:
            fprintf(filePointer, "boundary type: outflow\n"); 
            break;
        case SLIPWALL:
            fprintf(filePointer, "boundary type: slip wall\n"); 
            fprintf(filePointer, "temperature: %.6g\n", part->valueBC[n][ENTRYBC-1]);
            break;
        case NOSLIPWALL:
            fprintf(filePointer, "boundary type: noslip wall\n"); 
            fprintf(filePointer, "temperature: %.6g\n", part->valueBC[n][ENTRYBC-1]);
            break;
        case PERIODIC:
            fprintf(filePointer, "boundary type: periodic\n"); 
            break;
        default:
            FatalError("unidentified boundary type...");
            break;
    }
    return 0;
}
static int WriteInitializerData(FILE *filePointer, const Space *space, const int n)
{
    const Partition *part = &(space->part);
    switch (part->typeIC[n]) {
        case ICGLOBAL:
            break;
        case ICPLANE:
            fprintf(filePointer, "regional initialization: plane\n"); 
            fprintf(filePointer, "plane point x, y, z: %.6g, %.6g, %.6g\n", 
                    part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
            fprintf(filePointer, "plane normal nx, ny, nz: %.6g, %.6g, %.6g\n", 
                    part->valueIC[n][3], part->valueIC[n][4], part->valueIC[n][5]);
            break;
        case ICSPHERE:
            fprintf(filePointer, "regional initialization: sphere\n"); 
            fprintf(filePointer, "center point x, y, z: %.6g, %.6g, %.6g\n", 
                    part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
            fprintf(filePointer, "radius: %.6g\n", part->valueIC[n][6]);
            break;
        case ICBOX:
            fprintf(filePointer, "regional initialization: box\n"); 
            fprintf(filePointer, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n", 
                    part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
            fprintf(filePointer, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n", 
                    part->valueIC[n][3], part->valueIC[n][4], part->valueIC[n][5]);
            break;
        case ICCYLINDER:
            fprintf(filePointer, "regional initialization: cylinder\n"); 
            fprintf(filePointer, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n", 
                    part->valueIC[n][0], part->valueIC[n][1], part->valueIC[n][2]);
            fprintf(filePointer, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n", 
                    part->valueIC[n][3], part->valueIC[n][4], part->valueIC[n][5]);
            fprintf(filePointer, "radius: %.6g\n", part->valueIC[n][6]);
            break;
        default:
            break;
    }
    fprintf(filePointer, "density: %.6g\n", part->valueIC[n][ENTRYIC-5]);
    fprintf(filePointer, "x velocity: %.6g\n", part->valueIC[n][ENTRYIC-4]);
    fprintf(filePointer, "y velocity: %.6g\n", part->valueIC[n][ENTRYIC-3]);
    fprintf(filePointer, "z velocity: %.6g\n", part->valueIC[n][ENTRYIC-2]);
    fprintf(filePointer, "pressure: %.6g\n", part->valueIC[n][ENTRYIC-1]);
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
    fprintf(filePointer, "#                     Case Verification for ArtraCFD                          -\n");
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
    fprintf(filePointer, "termination time: %.6g\n", time->end); 
    fprintf(filePointer, "CFL condition number: %.6g\n", time->numCFL); 
    fprintf(filePointer, "maximum computing steps: %d\n", time->stepN); 
    fprintf(filePointer, "field data writing frequency: %d\n", time->writeN); 
    fprintf(filePointer, "data streamer: %d\n", time->dataStreamer); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Numerical Method <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "temporal scheme: %d\n", model->tScheme);
    fprintf(filePointer, "spatial scheme: %d\n", model->sScheme);
    fprintf(filePointer, "multidimensional method: %d\n", model->multidim);
    fprintf(filePointer, "Jacobian average: %d\n", model->jacobMean);
    fprintf(filePointer, "flux splitting method: %d\n", model->fluxSplit);
    fprintf(filePointer, "phase interaction: %d\n", model->fsi);
    fprintf(filePointer, "layers for reconstruction: %d\n", model->ibmLayer);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                       >> Material Properties <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "material: %d\n", model->mid); 
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
    fprintf(filePointer, "#                     >> Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    WriteInitializerData(filePointer, space, ICGLOBAL);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                     >> Boundary Condition <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian West\n"); 
    WriteBoundaryData(filePointer, space, PWB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian East\n"); 
    WriteBoundaryData(filePointer, space, PEB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian South\n"); 
    WriteBoundaryData(filePointer, space, PSB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian North\n"); 
    WriteBoundaryData(filePointer, space, PNB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Front\n"); 
    WriteBoundaryData(filePointer, space, PFB);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Back\n"); 
    WriteBoundaryData(filePointer, space, PBB);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                  >> Regional Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    for (int n = 1; n < part->countIC; ++n) {
        fprintf(filePointer, "#\n");
        WriteInitializerData(filePointer, space, n);
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Field Data Probes <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "point probe count: %d\n", time->pointProbeN);
    fprintf(filePointer, "line probe count: %d\n", time->lineProbeN);
    fprintf(filePointer, "curve probe count: %d\n", time->curveProbeN);
    fprintf(filePointer, "force probe count: %d\n", time->forceProbeN);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "point probe writing frequency: %d\n", time->pointWriteN);
    fprintf(filePointer, "line probe writing frequency: %d\n", time->lineWriteN);
    fprintf(filePointer, "body-conformal probe writing frequency: %d\n", time->curveWriteN);
    fprintf(filePointer, "surface force writing frequency: %d\n", time->forceWriteN);
    fprintf(filePointer, "#\n");
    for (int n = 0; n < time->pointProbeN; ++n) {
        fprintf(filePointer, "x, y, z of the point: %.6g, %.6g, %.6g\n", 
                time->pp[n][0], time->pp[n][1], time->pp[n][2]);
    }
    fprintf(filePointer, "#\n");
    for (int n = 0; n < time->lineProbeN; ++n) {
        fprintf(filePointer, "x, y, z of the first end point: %.6g, %.6g, %.6g\n", 
                time->lp[n][0], time->lp[n][1], time->lp[n][2]);
        fprintf(filePointer, "x, y, z of the second end point: %.6g, %.6g, %.6g\n", 
                time->lp[n][3], time->lp[n][4], time->lp[n][5]);
        fprintf(filePointer, "resolution: %.6g\n", time->lp[n][6]);
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
    const Real zero = 0.0;
    /* space */
    if ((zero >= (part->domain[X][MAX] - part->domain[X][MIN])) ||
            (zero >= (part->domain[Y][MAX] - part->domain[Y][MIN])) ||
            (zero >= (part->domain[Z][MAX] - part->domain[Z][MIN]))) {
        FatalError("wrong domian region values in case settings");
    }
    if ((1 > part->m[X]) || (1 > part->m[Y]) || (1 > part->m[Z])) {
        FatalError("too small mesh values in case settings");
    }
    /* time */
    if ((0 > time->restart) || (zero >= time->end) || (zero >= time->numCFL)) {
        FatalError("wrong values in time section of case settings");
    }
    /* numerical method */
    if ((0 > model->tScheme) || (0 > model->sScheme) || (0 > model->multidim) || 
            (0 > model->jacobMean) || (0 > model->fluxSplit) || (0 > model->fsi)) {
        FatalError("wrong values in numerical method of case settings");
    }
    /* material */
    if ((0 > model->mid)) {
        FatalError("wrong values in material section of case settings");
    }
    /* reference */
    if ((zero >= model->refL) || (zero >= model->refRho) || 
            (zero >= model->refV) || (zero >= model->refT)) {
        FatalError("wrong values in reference section of case settings");
    }
    return 0;
}
/* a good practice: end file with a newline */

