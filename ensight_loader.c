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
#include "ensight_stream.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "ensight.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int LoadEnsightCaseFile(EnsightSet *, Time *);
static int LoadEnsightVariableFile(Real *U, EnsightSet *,
        const Space *, const Partition *, const Flow *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Load necessary flow information from computed data
 */
int LoadComputedDataEnsight(Real *U, const Space *space, Time *time,
        const Partition *part, const Flow *flow)
{
    EnsightSet enSet = { /* initialize Ensight environment */
        .baseName = "ensight", /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .stringData = {'\0'}, /* string data recorder */
    };
    LoadEnsightCaseFile(&enSet, time);
    LoadEnsightVariableFile(U, &enSet, space, part, flow);
    return 0;
}
static int LoadEnsightCaseFile(EnsightSet *enSet, Time *time)
{
    FILE *filePointer = NULL;
    filePointer = fopen("restart.case", "r");
    if (NULL == filePointer) {
        FatalError("failed to open restart case file: restart.case...");
    }
    /* read information from file */
    char currentLine[200] = {'\0'}; /* store current line */
    /* set format specifier according to the type of Real */
    char format[25] = "%*s %*s %*s %*s %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(format, "%*s %*s %*s %*s %g", sizeof format); /* float type */
    }
    /* get rid of redundant lines */
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    /* get restart order number */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %*s %*s %d", &(time->outputCount)); 
    /* get restart time */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, format, &(time->currentTime)); 
    /* get current step number */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %*s %*s %d", &(time->stepCount)); 
    fclose(filePointer); /* close current opened file */
    /* store updated basename in filename */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s%05d", 
            enSet->baseName, time->outputCount); 
    /* basename is updated here! */
    snprintf(enSet->baseName, sizeof(EnsightString), "%s", enSet->fileName); 
    return 0;
}
static int LoadEnsightVariableFile(Real *U, EnsightSet *enSet,
        const Space *space, const Partition *part, const Flow *flow)
{
    FILE *filePointer = NULL;
    int idx = 0; /* linear array index math variable */
    EnsightReal data = 0.0; /* the ensight data format */
    const char nameSuffix[5][10] = {"rho", "u", "v", "w", "p"};
    for (int dim = 0; dim < DIMU; ++dim) {
        snprintf(enSet->fileName, sizeof(EnsightString), "%s.%s", enSet->baseName, nameSuffix[dim]);
        filePointer = fopen(enSet->fileName, "rb");
        if (NULL == filePointer) {
            FatalError("failed to open restart data files...");
        }
        fread(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        for (int partCount = 0, partNum = 1; partCount < NSUBPART; ++partCount) {
            fread(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            fread(&partNum, sizeof(int), 1, filePointer);
            fread(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            for (int k = part->kSub[partCount]; k < part->kSup[partCount]; ++k) {
                for (int j = part->jSub[partCount]; j < part->jSup[partCount]; ++j) {
                    for (int i = part->iSub[partCount]; i < part->iSup[partCount]; ++i) {
                        fread(&data, sizeof(EnsightReal), 1, filePointer);
                        idx = IndexMath(k, j, i, space) * DIMU;
                        switch (dim) {
                            case 0: /* rho */
                                U[idx] = data;
                                break;
                            case 1: /* u */
                                U[idx+1] = U[idx] * data;
                                break;
                            case 2: /* v */
                                U[idx+2] = U[idx] * data;
                                break;
                            case 3: /* w */
                                U[idx+3] = U[idx] * data;
                                break;
                            case 4: /* p */
                                U[idx+4] = 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx] + data / (flow->gamma - 1.0);
                                break;
                            default:
                                break;
                        }
                    }
                }
            }
        }
        fclose(filePointer); /* close current opened file */
    }
    return 0;
}
/* a good practice: end file with a newline */

