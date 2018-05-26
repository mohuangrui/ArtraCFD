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
#include "ensight.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ReadCaseFile(Time *, EnsightSet *);
static int ReadStructuredData(Space *, const Model *, EnsightSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ReadStructuredDataEnsight(Time *time, Space *space, const Model *model)
{
    EnsightSet enSet = { /* initialize environment */
        .rootName = "field", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .stringData = {'\0'}, /* string data recorder */
    };
    snprintf(enSet.baseName, sizeof(EnsightString), "%s%05d", 
            enSet.rootName, time->writeC); 
    ReadCaseFile(time, &enSet);
    ReadStructuredData(space, model, &enSet);
    return 0;
}
static int ReadCaseFile(Time *time, EnsightSet *enSet)
{
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.case", 
            enSet->baseName); 
    FILE *filePointer = fopen(enSet->fileName, "r");
    if (NULL == filePointer) {
        FatalError("failed to open case file...");
    }
    /* read information from file */
    String currentLine = {'\0'}; /* store current line */
    int nscan = 0; /* read conversion count */
    ReadInLine(filePointer, "VARIABLE");
    /* set format specifier according to the type of Real */
    char format[25] = "%*s %*s %*s %*s %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(format, "%*s %*s %*s %*s %g", sizeof format); /* float type */
    }
    Fgets(currentLine, sizeof currentLine, filePointer);
    nscan = sscanf(currentLine, format, &(time->now)); 
    VerifyReadConversion(nscan, 1);
    Fgets(currentLine, sizeof currentLine, filePointer);
    nscan = sscanf(currentLine, "%*s %*s %*s %*s %d", &(time->stepC)); 
    VerifyReadConversion(nscan, 1);
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int ReadStructuredData(Space *space, const Model *model, EnsightSet *enSet)
{
    FILE *filePointer = NULL;
    EnsightReal data = 0.0; /* the Ensight data format */
    const char scalar[5][5] = {"rho", "u", "v", "w", "p"};
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    for (int count = 0; count < DIMU; ++count) {
        snprintf(enSet->fileName, sizeof(EnsightString), "%s.%s", enSet->baseName, scalar[count]);
        filePointer = fopen(enSet->fileName, "rb");
        if (NULL == filePointer) {
            FatalError("failed to open data file...");
        }
        Fread(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        for (int p = PIN, partNum = 1; p < NPARTWRITE; ++p) {
            Fread(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            Fread(&partNum, sizeof(int), 1, filePointer);
            Fread(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
                for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
                    for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                        idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                        U = node[idx].U[TO];
                        Fread(&data, sizeof(EnsightReal), 1, filePointer);
                        switch (count) {
                            case 0: /* rho */
                                U[0] = data;
                                break;
                            case 1: /* u */
                                U[1] = U[0] * data;
                                break;
                            case 2: /* v */
                                U[2] = U[0] * data;
                                break;
                            case 3: /* w */
                                U[3] = U[0] * data;
                                break;
                            case 4: /* p */
                                U[4] = 0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / U[0] + 
                                    data / (model->gamma - 1.0);
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

