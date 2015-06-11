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
#include "parasight_stream.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "ensight.h" 
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int InitializeParasightTransientCaseFile(EnsightSet *);
static int WriteParasightCaseFile(const Time *, EnsightSet *);
static int WriteParasightGeometryFile(const Space *, const Partition *, EnsightSet *);
static int WriteParasightVariableFile(const Real *, const Space *, const Model *,
        const Partition *, EnsightSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function write computed data to files with Parasight data format, 
 * including transient and steady output with file names consists of the 
 * default base file name and export step tag. 
 */
int WriteComputedDataParasight(const Real *U, const Space *space, 
        const Time *time, const Model *model, const Partition *part)
{
    ShowInformation("  writing field data to file...");
    EnsightSet enSet = { /* initialize Parasight environment */
        .baseName = "parasight", /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .stringData = {'\0'}, /* string data recorder */
    };
    if (0 == time->stepCount) { /* this is the initialization step */
        InitializeParasightTransientCaseFile(&enSet);
        WriteParasightGeometryFile(space, part, &enSet);
    }
    WriteParasightCaseFile(time, &enSet);
    WriteParasightVariableFile(U, space, model, part, &enSet);
    return 0;
}
/*
 * Parasight transient case file
 * This function initializes an overall transient case file.
 */
int InitializeParasightTransientCaseFile(EnsightSet *enSet)
{
    FILE *filePointer = fopen("parasight.case", "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to transient case file...");
    }
    /* output information to file */
    fprintf(filePointer, "FORMAT\n"); 
    fprintf(filePointer, "type: ensight gold\n"); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "GEOMETRY\n"); 
    fprintf(filePointer, "model:  parasight.geo\n"); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "VARIABLE\n"); 
    fprintf(filePointer, "scalar per node:  1  rho  %s*****.rho\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:  1  u    %s*****.u\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:  1  v    %s*****.v\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:  1  w    %s*****.w\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:  1  p    %s*****.p\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:  1  T    %s*****.T\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:  1  id   %s*****.id\n", enSet->baseName); 
    fprintf(filePointer, "vector per node:  1  Vel  %s*****.Vel\n", enSet->baseName); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "TIME\n"); 
    fprintf(filePointer, "time set:         1\n"); 
    fprintf(filePointer, "number of steps:          0          \n");
    fprintf(filePointer, "filename start number:    0\n"); 
    fprintf(filePointer, "filename increment:       1\n"); 
    fprintf(filePointer, "time values:  "); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * Parasight case files
 */
static int WriteParasightCaseFile(const Time *time, EnsightSet *enSet)
{
    /*
     * Write the steady case file of current step.
     * To get the target file name, the function snprintf is used. It prints 
     * specified format to a string (char array) with buffer overflow 
     * protection. 
     * NOTE: if memeory locations of input objects overlap, the behavior of
     * snprintf is undefined!
     * To make all step number has the same digit number, in the
     * format specifiers, zero padding is added to the integer value. a zero
     * '0' character indicating that zero-padding should be used rather than
     * blank-padding.
     */
    /* store updated basename in filename */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s%05d", 
            enSet->baseName, time->outputCount); 
    /* basename is updated here! */
    snprintf(enSet->baseName, sizeof(EnsightString), "%s", enSet->fileName); 
    /* current filename */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.case", enSet->baseName); 
    FILE *filePointer = fopen(enSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to steady case file...");
    }
    /* output information to file */
    fprintf(filePointer, "FORMAT\n"); 
    fprintf(filePointer, "type: ensight gold\n"); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "GEOMETRY\n"); 
    fprintf(filePointer, "model:  parasight.geo\n"); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "VARIABLE\n"); 
    fprintf(filePointer, "constant per case:  Order %d\n", time->outputCount);
    fprintf(filePointer, "constant per case:  Time  %.6g\n", time->now);
    fprintf(filePointer, "constant per case:  Step  %d\n", time->stepCount);
    fprintf(filePointer, "scalar per node:    rho   %s.rho\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    u     %s.u\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    v     %s.v\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    w     %s.w\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    p     %s.p\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    T     %s.T\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    id    %s.id\n", enSet->baseName); 
    fprintf(filePointer, "vector per node:    Vel   %s.Vel\n", enSet->baseName); 
    fprintf(filePointer, "\n"); 
    fclose(filePointer); /* close current opened file */
    /*
     * Add information to the transient case file
     */
    /* correct the number of steps in transient case */
    filePointer = fopen("parasight.case", "r+");
    if (NULL == filePointer) {
        FatalError("failed to add data to transient file...");
    }
    /* seek the target line for adding information */
    char currentLine[200] = {'\0'}; /* store the current read line */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "time set", 8)) {
            break;
        }
    }
    fprintf(filePointer, "number of steps:          %d", (time->outputCount + 1)); 
    /* add the time flag of current export to the transient case */
    fseek(filePointer, 0, SEEK_END); // seek to the end of file
    if ((time->outputCount % 5) == 0) { /* print to a new line every x outputs */
        fprintf(filePointer, "\n"); 
    }
    fprintf(filePointer, "%.6g ", time->now); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * Parasight geometry file
 */
static int WriteParasightGeometryFile(const Space *space, const Partition *part, EnsightSet *enSet)
{
    /*
     * Write the geometry file (Binary Form).
     * Parasight Maximums: maximum number of nodes in a part is 2GB.
     */
    snprintf(enSet->fileName, sizeof(EnsightString), "parasight.geo");
    FILE *filePointer = fopen(enSet->fileName, "wb");
    if (NULL == filePointer) {
        FatalError("failed to write geometry file...");
    }
    /*
     * Output information to file, need to strictly follow the Parasight data format.
     * NOTE: if memeory locations of input objects overlap, the behavior of
     * strncpy is undefined!
     * In fwrite, the first size is the sizeof an object, which is given in the
     * units of chars, And the second size (count) is the number of object 
     * that need to be written.
     */
    int nodeCount[3] = {0, 0, 0}; /* i j k node number in part */
    const int partNum = 1; /* only one part is needed */
    EnsightReal data = 0.0; /* the Parasight data format */
    /* description at the beginning */
    strncpy(enSet->stringData, "C Binary", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "Parasight Geometry File", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "Written by ArtraCFD", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    /* node id and extents settings */
    strncpy(enSet->stringData, "node id off", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "element id off", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "part", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    fwrite(&partNum, sizeof(int), 1, filePointer);
    strncpy(enSet->stringData, "entire domain", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "block", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    /* this line is the total number of nodes in i, j, k */
    nodeCount[0] = (part->iSup[0] - part->iSub[0]); 
    nodeCount[1] = (part->jSup[0] - part->jSub[0]);
    nodeCount[2] = (part->kSup[0] - part->kSub[0]);
    fwrite(nodeCount, sizeof(int), 3, filePointer);
    /* now output the x coordinates of all nodes in current part */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                data = ComputeX(i, space);
                fwrite(&data, sizeof(EnsightReal), 1, filePointer);
            }
        }
    }
    /* now output the y coordinates of all nodes in current part */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                data = ComputeY(j, space);
                fwrite(&data, sizeof(EnsightReal), 1, filePointer);
            }
        }
    }
    /* now output the z coordinates of all nodes in current part */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                data = ComputeZ(k, space);
                fwrite(&data, sizeof(EnsightReal), 1, filePointer);
            }
        }
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * Parasight variables files
 * The values for each node of the structured block are output in 
 * the same IJK order as the coordinates. (The number of nodes in the
 * part are obtained from the corresponding Parasight Gold geometry file.)
 */
static int WriteParasightVariableFile(const Real *U, const Space *space,
        const Model *model, const Partition *part, EnsightSet *enSet)
{
    FILE *filePointer = NULL;
    int idx = 0; /* linear array index math variable */
    EnsightReal data = 0.0; /* the Parasight data format */
    /*
     * Write the scalar field (Binary Form)
     */
    const char nameSuffix[7][5] = {"rho", "u", "v", "w", "p", "T", "id"};
    const int partNum = 1;
    for (int dim = 0; dim < 7; ++dim) {
        snprintf(enSet->fileName, sizeof(EnsightString), "%s.%s", enSet->baseName, nameSuffix[dim]);
        filePointer = fopen(enSet->fileName, "wb");
        if (NULL == filePointer) {
            FatalError("failed to write data file...");
        }
        /* first line description per file */
        strncpy(enSet->stringData, "scalar variable", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        /* binary file format */
        strncpy(enSet->stringData, "part", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        fwrite(&partNum, sizeof(int), 1, filePointer);
        strncpy(enSet->stringData, "block", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        /* now output the scalar value at each node in current part */
        for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
            for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
                for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                    idx = IndexMath(k, j, i, space) * DIMU;
                    switch (dim) {
                        case 0: /* rho */
                            data = U[idx];
                            break;
                        case 1: /* u */
                            data = U[idx+1] / U[idx];
                            break;
                        case 2: /* v */
                            data = U[idx+2] / U[idx];
                            break;
                        case 3: /* w */
                            data = U[idx+3] / U[idx];
                            break;
                        case 4: /* p */
                            data = ComputePressure(idx, U, model);
                            break;
                        case 5: /* T */
                            data = ComputeTemperature(idx, U, model);
                            break;
                        case 6: /* node flag */
                            idx = idx / DIMU;
                            data = space->nodeFlag[idx];
                        default:
                            break;
                    }
                    fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                }
            }
        }
        fclose(filePointer); /* close current opened file */
    }
    /*
     * Write the velocity vector field (Binary Form)
     */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.Vel", enSet->baseName);
    filePointer = fopen(enSet->fileName, "wb");
    if (NULL == filePointer) {
        FatalError("failed to write data file...");
    }
    /* binary file format */
    strncpy(enSet->stringData, "vector variable", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "part", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    fwrite(&partNum, sizeof(int), 1, filePointer);
    strncpy(enSet->stringData, "block", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    /*
     * Now output the vector components at each node in current part
     * dimension index of u, v, w is 1, 2, 3 in U in each part, 
     * write u, v, w sequentially
     */
    for (int dim = 1; dim < 4; ++dim) {
        for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
            for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
                for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                    idx = IndexMath(k, j, i, space) * DIMU;
                    switch (dim) {
                        case 1: /* u */
                            data = U[idx+1] / U[idx];
                            break;
                        case 2: /* v */
                            data = U[idx+2] / U[idx];
                            break;
                        case 3: /* w */
                            data = U[idx+3] / U[idx];
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
    return 0;
}
/* a good practice: end file with a newline */

