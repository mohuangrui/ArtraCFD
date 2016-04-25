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
#include "paraview.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ReadCaseFile(Time *, ParaviewSet *);
static int ReadStructuredData(Space *, const Model *, ParaviewSet *);
static int PointPolyDataReader(Geometry *, const Time *);
static int ReadPointPolyData(Geometry *, ParaviewSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ReadStructuredDataParaview(Space *space, Time *time, const Model *model)
{
    ParaviewSet paraSet = { /* initialize ParaviewSet environment */
        .rootName = "field", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .fileExt = ".vts", /* data file extension */
        .intType = "Int32", /* paraview int type */
        .floatType = "Float32", /* paraview float type */
        .byteOrder = "LittleEndian" /* byte order of data */
    };
    snprintf(paraSet.baseName, sizeof(ParaviewString), "%s%05d", 
            paraSet.rootName, time->countOutput); 
    ReadCaseFile(time, &paraSet);
    ReadStructuredData(space, model, &paraSet);
    return 0;
}
static int ReadCaseFile(Time *time, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.pvd", 
            paraSet->baseName); 
    FILE *filePointer = fopen(paraSet->fileName, "r");
    if (NULL == filePointer) {
        FatalError("failed to open case file...");
    }
    /* read information from file */
    String currentLine = {'\0'}; /* store current line */
    /* get rid of redundant lines */
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    /* get restart time */
    /* set format specifier according to the type of Real */
    char format[25] = {'\0'}; /* format information */
    strncpy(format, "%*s %*s %lg", sizeof format); /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(format, "%*s %*s %g", sizeof format); /* float type */
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, format, &(time->now)); 
    /* get current step number */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %d", &(time->countStep)); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int ReadStructuredData(Space *space, const Model *model, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%s", paraSet->baseName, paraSet->fileExt); 
    FILE *filePointer = fopen(paraSet->fileName, "r");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    ParaviewReal data = 0.0; /* paraview scalar data */
    /* set format specifier according to the type of Real */
    char format[5] = "%lg"; /* default is double type */
    if (sizeof(ParaviewReal) == sizeof(float)) {
        strncpy(format, "%g", sizeof format); /* float type */
    }
    int idx = 0; /* linear array index math variable */
    Node *node = space->node;
    Real *restrict U = NULL;
    const Partition *restrict part = &(space->part);
    /* get rid of redundant lines */
    String currentLine = {'\0'}; /* store current line */
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    for (int count = 0; count < DIMU; ++count) {
        fgets(currentLine, sizeof currentLine, filePointer);
        for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
            for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
                for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    U = node[idx].U[C];
                    fscanf(filePointer, format, &data);
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
                            U[4] = 0.5 * (U[1] * U[1] + U[2] * U[2] + U[3] * U[3]) / U[0] + data / (model->gamma - 1.0);
                            break;
                        default:
                            break;
                    }
                }
            }
        }
        fgets(currentLine, sizeof currentLine, filePointer); /* get rid of the end of line of data */
        fgets(currentLine, sizeof currentLine, filePointer);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int PointPolyDataReader(Geometry *geo, const Time *time)
{
    ParaviewSet paraSet = { /* initialize ParaviewSet environment */
        .rootName = "geo_sph", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .fileExt = ".vtp", /* data file extension */
        .intType = "Int32", /* paraview int type */
        .floatType = "Float32", /* paraview float type */
        .byteOrder = "LittleEndian" /* byte order of data */
    };
    snprintf(paraSet.baseName, sizeof(ParaviewString), "%s%05d", 
            paraSet.rootName, time->countOutput); 
    ReadPointPolyData(geo, &paraSet);
    return 0;
}
static int ReadPointPolyData(Geometry *geo, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%s", paraSet->baseName, paraSet->fileExt); 
    FILE *filePointer = fopen(paraSet->fileName, "r");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    ParaviewReal data = 0.0; /* paraview scalar data */
    /* set format specifier according to the type of Real */
    char format[5] = "%lg"; /* default is double type */
    if (sizeof(ParaviewReal) == sizeof(float)) {
        strncpy(format, "%g", sizeof format); /* float type */
    }
    /* get rid of redundant lines */
    String currentLine = {'\0'}; /* store current line */
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
}
/* a good practice: end file with a newline */

