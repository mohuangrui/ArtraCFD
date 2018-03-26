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
#include <float.h> /* size of floating point values */
#include "computational_geometry.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ReadCaseFile(Time *, ParaviewSet *);
static int ReadStructuredData(Space *, const Model *, ParaviewSet *);
static int PointPolyDataReader(const Time *, Geometry *);
static int ReadPointPolyData(const int, const int, Geometry *, ParaviewSet *);
static int PolygonPolyDataReader(const Time *, Geometry *);
static int ReadPolygonPolyData(const int, const int, Geometry *, ParaviewSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ReadStructuredDataParaview(Time *time, Space *space, const Model *model)
{
    ParaviewSet paraSet = { /* initialize environment */
        .rootName = "field", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .fileExt = ".vts", /* data file extension */
        .intType = "Int32", /* paraview int type */
        .floatType = "Float32", /* paraview float type */
        .byteOrder = "LittleEndian" /* byte order of data */
    };
    snprintf(paraSet.baseName, sizeof(ParaviewString), "%s%05d", 
            paraSet.rootName, time->writeC); 
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
    int nscan = 0; /* read conversion count */
    ReadInLine(filePointer, "<!--");
    /* set format specifier according to the type of Real */
    char format[10] = {'\0'}; /* format information */
    strncpy(format, "%*s %lg", sizeof format); /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(format, "%*s %g", sizeof format); /* float type */
    }
    Fgets(currentLine, sizeof currentLine, filePointer);
    nscan = sscanf(currentLine, format, &(time->now)); 
    VerifyReadConversion(nscan, 1);
    Fgets(currentLine, sizeof currentLine, filePointer);
    nscan = sscanf(currentLine, "%*s %d", &(time->stepC)); 
    VerifyReadConversion(nscan, 1);
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
    int nscan = 0; /* read conversion count */
    ParaviewReal data = 0.0; /* paraview scalar data */
    /* set format specifier according to the type of Real */
    char format[5] = "%lg"; /* default is double type */
    if (sizeof(ParaviewReal) == sizeof(float)) {
        strncpy(format, "%g", sizeof format); /* float type */
    }
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    /* get rid of redundant lines */
    String currentLine = {'\0'}; /* store current line */
    ReadInLine(filePointer, "<PointData>");
    for (int count = 0; count < DIMU; ++count) {
        Fgets(currentLine, sizeof currentLine, filePointer);
        for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
            for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
                for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    U = node[idx].U[TO];
                    nscan = fscanf(filePointer, format, &data);
                    VerifyReadConversion(nscan, 1);
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
        Fgets(currentLine, sizeof currentLine, filePointer); /* get rid of the end of line of data */
        Fgets(currentLine, sizeof currentLine, filePointer);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
int ReadPolyDataParaview(const Time *time, Geometry *geo)
{
    if (0 != geo->sphN) {
        PointPolyDataReader(time, geo);
    }
    if (0 != geo->stlN) {
        PolygonPolyDataReader(time, geo);
    }
    return 0;
}
static int PointPolyDataReader(const Time *time, Geometry *geo)
{
    ParaviewSet paraSet = { /* initialize environment */
        .rootName = "geo_sph", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .fileExt = ".vtp", /* data file extension */
        .intType = "Int32", /* paraview int type */
        .floatType = "Float32", /* paraview float type */
        .byteOrder = "LittleEndian" /* byte order of data */
    };
    snprintf(paraSet.baseName, sizeof(ParaviewString), "%s%05d", 
            paraSet.rootName, time->writeC); 
    ReadPointPolyData(0, geo->sphN, geo, &paraSet);
    return 0;
}
static int ReadPointPolyData(const int start, const int end, Geometry *geo, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%s", paraSet->baseName, paraSet->fileExt); 
    FILE *filePointer = fopen(paraSet->fileName, "r");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    ReadInLine(filePointer, "<!--");
    ReadPolyhedronStateData(start, end, filePointer, geo);
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int PolygonPolyDataReader(const Time *time, Geometry *geo)
{
    ParaviewSet paraSet = { /* initialize environment */
        .rootName = "geo_stl", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .fileExt = ".vtp", /* data file extension */
        .intType = "Int32", /* paraview int type */
        .floatType = "Float32", /* paraview float type */
        .byteOrder = "LittleEndian" /* byte order of data */
    };
    snprintf(paraSet.baseName, sizeof(ParaviewString), "%s%05d", 
            paraSet.rootName, time->writeC); 
    ReadPolygonPolyData(geo->sphN, geo->totN, geo, &paraSet);
    return 0;
}
static int ReadPolygonPolyData(const int start, const int end, Geometry *geo, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%s", paraSet->baseName, paraSet->fileExt); 
    FILE *filePointer = fopen(paraSet->fileName, "r");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    int nscan = 0; /* read conversion count */
    ParaviewReal Vec[3] = {0.0}; /* paraview vector data */
    Polyhedron *poly  = NULL;
    /* set format specifier according to the type of Real */
    char format[15] = "%lg %lg %lg"; /* default is double type */
    if (sizeof(ParaviewReal) == sizeof(float)) {
        strncpy(format, "%g %g %g", sizeof format); /* float type */
    }
    /* get rid of redundant lines */
    String currentLine = {'\0'}; /* store current line */
    ReadInLine(filePointer, "<PolyData>");
    for (int m = start; m < end; ++m) {
        poly = geo->poly + m;
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        nscan = sscanf(currentLine, "%*s %*s %d", &(poly->vertN)); 
        VerifyReadConversion(nscan, 1);
        Fgets(currentLine, sizeof currentLine, filePointer);
        nscan = sscanf(currentLine, "%*s %*s %d", &(poly->edgeN)); 
        VerifyReadConversion(nscan, 1);
        Fgets(currentLine, sizeof currentLine, filePointer);
        nscan = sscanf(currentLine, "%*s %*s %d", &(poly->faceN)); 
        VerifyReadConversion(nscan, 1);
        Fgets(currentLine, sizeof currentLine, filePointer);
        AllocatePolyhedronMemory(poly->vertN, poly->edgeN, poly->faceN, poly);
        poly->edgeN = 0; /* reset edge count before applying edge adding */
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        for (int n = 0; n < poly->vertN; ++n) {
            nscan = fscanf(filePointer, format, &(Vec[X]), &(Vec[Y]), &(Vec[Z]));
            VerifyReadConversion(nscan, 3);
            poly->v[n][X] = Vec[X];
            poly->v[n][Y] = Vec[Y];
            poly->v[n][Z] = Vec[Z];
        }
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        Fgets(currentLine, sizeof currentLine, filePointer);
        for (int n = 0; n < poly->faceN; ++n) {
            nscan = fscanf(filePointer, "%d %d %d", &(poly->f[n][0]), 
                    &(poly->f[n][1]), &(poly->f[n][2]));
            VerifyReadConversion(nscan, 3);
            AddEdge(poly->f[n][0], poly->f[n][1], n, poly); 
            AddEdge(poly->f[n][1], poly->f[n][2], n, poly); 
            AddEdge(poly->f[n][2], poly->f[n][0], n, poly); 
        }
        QuickSortEdge(poly->edgeN, poly->e);
        ReadInLine(filePointer, "</Piece>");
    }
    ReadInLine(filePointer, "<!--");
    ReadPolyhedronStateData(start, end, filePointer, geo);
    fclose(filePointer); /* close current opened file */
    return 0;
}
int ReadPolyhedronStateData(const int start, const int end, FILE *filePointer, Geometry *geo)
{
    String currentLine = {'\0'}; /* store the current read line */
    int nscan = 0; /* read conversion count */
    /* set format specifier according to the type of Real */
    char formatI[100] = "%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %d";
    char formatJ[100] = "%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg";
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %d", sizeof formatI);
        strncpy(formatJ, "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g", sizeof formatJ);
    }
    Polyhedron *poly  = NULL;
    const Real zero = 0.0;
    for (int n = start; n < end; ++n) {
        poly = geo->poly + n;
        Fgets(currentLine, sizeof currentLine, filePointer);
        nscan = sscanf(currentLine, formatI,
                &(poly->O[X]), &(poly->O[Y]), &(poly->O[Z]), &(poly->r),
                &(poly->V[TO][X]), &(poly->V[TO][Y]), &(poly->V[TO][Z]),
                &(poly->W[TO][X]), &(poly->W[TO][Y]), &(poly->W[TO][Z]),
                &(poly->rho), &(poly->T), &(poly->cf),
                &(poly->area), &(poly->volume), &(poly->mid));
        VerifyReadConversion(nscan, 16);
        Fgets(currentLine, sizeof currentLine, filePointer);
        nscan = sscanf(currentLine, formatJ,
                &(poly->at[TO][X]), &(poly->at[TO][Y]), &(poly->at[TO][Z]),
                &(poly->ar[TO][X]), &(poly->ar[TO][Y]), &(poly->ar[TO][Z]),
                &(poly->at[TN][X]), &(poly->at[TN][Y]), &(poly->at[TN][Z]),
                &(poly->g[X]), &(poly->g[Y]), &(poly->g[Z]),
                &(poly->ar[TN][X]), &(poly->ar[TN][Y]), &(poly->ar[TN][Z]),
                &(poly->to));
        VerifyReadConversion(nscan, 16);
        if (zero >= poly->to) {
            poly->to = FLT_MAX;
        }
        if (geo->sphN > n) {
            poly->faceN = 0; /* analytical sphere tag */
            poly->facet = NULL;
        }
        for (int s = 0; s < DIMS; ++s) {
            poly->V[TN][s] = poly->V[TO][s];
            poly->W[TN][s] = poly->W[TO][s];
        }
    }
    return 0;
}
/* a good practice: end file with a newline */

