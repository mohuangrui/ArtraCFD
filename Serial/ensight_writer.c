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
static int InitializeTransientCaseFile(EnsightSet *);
static int WriteCaseFile(const Time *, EnsightSet *);
static int WriteGeometryFile(const Space *, EnsightSet *);
static int WriteStructuredData(const Space *, const Model *, EnsightSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteStructuredDataEnsight(const Time *time, const Space *space, const Model *model)
{
    EnsightSet enSet = { /* initialize environment */
        .rootName = "field", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .stringData = {'\0'}, /* string data recorder */
    };
    snprintf(enSet.baseName, sizeof(EnsightString), "%s%05d", 
            enSet.rootName, time->writeC); 
    if (0 == time->stepC) { /* this is the initialization step */
        InitializeTransientCaseFile(&enSet);
        WriteGeometryFile(space, &enSet);
    }
    WriteCaseFile(time, &enSet);
    WriteStructuredData(space, model, &enSet);
    return 0;
}
static int InitializeTransientCaseFile(EnsightSet *enSet)
{
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.case", 
            enSet->rootName); 
    FILE *filePointer = fopen(enSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to initialize transient case file...");
    }
    /* output information to file */
    fprintf(filePointer, "FORMAT\n"); 
    fprintf(filePointer, "type: ensight gold\n"); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "GEOMETRY\n"); 
    fprintf(filePointer, "model:  %s.geo\n", enSet->rootName); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "VARIABLE\n"); 
    fprintf(filePointer, "scalar per node:  1  rho  %s*****.rho\n", enSet->rootName); 
    fprintf(filePointer, "scalar per node:  1  u    %s*****.u\n", enSet->rootName); 
    fprintf(filePointer, "scalar per node:  1  v    %s*****.v\n", enSet->rootName); 
    fprintf(filePointer, "scalar per node:  1  w    %s*****.w\n", enSet->rootName); 
    fprintf(filePointer, "scalar per node:  1  p    %s*****.p\n", enSet->rootName); 
    fprintf(filePointer, "scalar per node:  1  T    %s*****.T\n", enSet->rootName); 
    fprintf(filePointer, "scalar per node:  1  gid  %s*****.gid\n", enSet->rootName); 
    fprintf(filePointer, "vector per node:  1  Vel  %s*****.Vel\n", enSet->rootName); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "TIME\n"); 
    fprintf(filePointer, "time set: 1\n"); 
    fprintf(filePointer, "number of steps:          0          \n");
    fprintf(filePointer, "filename start number:    0\n"); 
    fprintf(filePointer, "filename increment:       1\n"); 
    fprintf(filePointer, "time values:  "); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteCaseFile(const Time *time, EnsightSet *enSet)
{
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.case", 
            enSet->baseName); 
    FILE *filePointer = fopen(enSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to open case file...");
    }
    /* output information to file */
    fprintf(filePointer, "FORMAT\n"); 
    fprintf(filePointer, "type: ensight gold\n"); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "GEOMETRY\n"); 
    fprintf(filePointer, "model:  %s.geo\n", enSet->rootName); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "VARIABLE\n"); 
    fprintf(filePointer, "constant per case:  Time  %.6g\n", time->now);
    fprintf(filePointer, "constant per case:  Step  %d\n", time->stepC);
    fprintf(filePointer, "scalar per node:    rho   %s.rho\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    u     %s.u\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    v     %s.v\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    w     %s.w\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    p     %s.p\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    T     %s.T\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    gid   %s.gid\n", enSet->baseName); 
    fprintf(filePointer, "vector per node:    Vel   %s.Vel\n", enSet->baseName); 
    fprintf(filePointer, "\n"); 
    fclose(filePointer); /* close current opened file */
    /*
     * Add information to the transient case file
     */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.case", 
            enSet->rootName); 
    filePointer = fopen(enSet->fileName, "r+");
    if (NULL == filePointer) {
        FatalError("failed to add data to transient file...");
    }
    /* seek the target line for adding information */
    ReadInLine(filePointer, "time set: 1");
    fprintf(filePointer, "number of steps:          %d", (time->writeC + 1)); 
    /* add the time flag of current export to the transient case */
    fseek(filePointer, 0, SEEK_END); // seek to the end of file
    if ((time->writeC % 5) == 0) { /* print to a new line every x outputs */
        fprintf(filePointer, "\n"); 
    }
    fprintf(filePointer, "%.6g ", time->now); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteGeometryFile(const Space *space, EnsightSet *enSet)
{
    /*
     * Write the geometry file in Binary Form.
     * Maximums: maximum number of nodes in a part is 2GB.
     */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.geo", 
            enSet->rootName); 
    FILE *filePointer = fopen(enSet->fileName, "wb");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    EnsightReal data = 0.0; /* the Ensight data format */
    const Partition *restrict part = &(space->part);
    IntVec nodeCount = {0}; /* i, j, k node number in each part */
    /*
     * Output information to file, need to strictly follow the Ensight data format.
     * In fwrite, the first size is the sizeof an object, which is given in the
     * units of chars, And the second size (count) is the number of object 
     * that need to be written.
     */
    /* description at the beginning */
    strncpy(enSet->stringData, "C Binary", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "Ensight Geometry File", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "Written by ArtraCFD", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    /* node id and extents settings */
    strncpy(enSet->stringData, "node id off", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "element id off", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    for (int p = PIN, partNum = 1; p < NPARTWRITE; ++p, ++partNum) {
        strncpy(enSet->stringData, "part", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        fwrite(&partNum, sizeof(int), 1, filePointer);
        snprintf(enSet->stringData, sizeof(EnsightString), "part %d", p);
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        strncpy(enSet->stringData, "block", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        nodeCount[X] = part->ns[p][X][MAX] - part->ns[p][X][MIN]; 
        nodeCount[Y] = part->ns[p][Y][MAX] - part->ns[p][Y][MIN]; 
        nodeCount[Z] = part->ns[p][Z][MAX] - part->ns[p][Z][MIN]; 
        fwrite(nodeCount, sizeof(int), 3, filePointer);
        /* now output the x coordinates of all nodes in current part */
        for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
            for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
                for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                    data = PointSpace(i, part->domain[X][MIN], part->d[X], part->ng);
                    fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                }
            }
        }
        /* now output the y coordinates of all nodes in current part */
        for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
            for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
                for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                    data = PointSpace(j, part->domain[Y][MIN], part->d[Y], part->ng);
                    fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                }
            }
        }
        /* now output the z coordinates of all nodes in current part */
        for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
            for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
                for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                    data = PointSpace(k, part->domain[Z][MIN], part->d[Z], part->ng);
                    fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                }
            }
        }
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * The values for each node of the structured block are output in 
 * the same IJK order as the coordinates. (The number of nodes in the
 * part are obtained from the corresponding geometry file.)
 */
static int WriteStructuredData(const Space *space, const Model *model, EnsightSet *enSet)
{
    FILE *filePointer = NULL;
    EnsightReal data = 0.0; /* the Ensight data format */
    const char scalar[7][5] = {"rho", "u", "v", "w", "p", "T", "gid"};
    const Partition *restrict part = &(space->part);
    const Node *const node = space->node;
    const Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    for (int count = 0; count < 7; ++count) {
        snprintf(enSet->fileName, sizeof(EnsightString), "%s.%s", enSet->baseName, scalar[count]);
        filePointer = fopen(enSet->fileName, "wb");
        if (NULL == filePointer) {
            FatalError("failed to open data file...");
        }
        /* first line description per file */
        strncpy(enSet->stringData, "scalar variable", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        for (int p = PIN, partNum = 1; p < NPARTWRITE; ++p, ++partNum) {
            /* binary file format */
            strncpy(enSet->stringData, "part", sizeof(EnsightString));
            fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            fwrite(&partNum, sizeof(int), 1, filePointer);
            strncpy(enSet->stringData, "block", sizeof(EnsightString));
            fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            /* now output the scalar value at each node in current part */
            for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
                for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
                    for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                        idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                        U = node[idx].U[TO];
                        switch (count) {
                            case 0: /* rho */
                                data = U[0];
                                break;
                            case 1: /* u */
                                data = U[1] / U[0];
                                break;
                            case 2: /* v */
                                data = U[2] / U[0];
                                break;
                            case 3: /* w */
                                data = U[3] / U[0];
                                break;
                            case 4: /* p */
                                data = ComputePressure(model->gamma, U);
                                break;
                            case 5: /* T */
                                data = ComputeTemperature(model->cv, U);
                                break;
                            case 6: /* node flag */
                                data = node[idx].gid;
                                break;
                            default:
                                break;
                        }
                        fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                    }
                }
            }
        }
        fclose(filePointer); /* close current opened file */
    }
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.Vel", enSet->baseName);
    filePointer = fopen(enSet->fileName, "wb");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    /* binary file format */
    strncpy(enSet->stringData, "vector variable", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    for (int p = PIN, partNum = 1; p < NPARTWRITE; ++p, ++partNum) {
        strncpy(enSet->stringData, "part", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        fwrite(&partNum, sizeof(int), 1, filePointer);
        strncpy(enSet->stringData, "block", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        for (int count = 1; count < 4; ++count) {
            for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
                for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
                    for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                        idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                        U = node[idx].U[TO];
                        data = U[count] / U[0];
                        fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                    }
                }
            }
        }
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

