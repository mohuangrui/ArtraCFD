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
#include "geometry_stream.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "stl.h"
#include "paraview.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NonrestartGeometryReader(Geometry *);
static int RestartGeometryReader(const Time *, Geometry *);
static int ReadSphereFile(const char *, Geometry *);
static int ReadPolygonStatusData(FILE **, Polygon *);
static int WritePolygonStatusData(FILE **, Polygon *);
static int ReadGeometryDataParaview(const Time *, Geometry *);
static int WriteGeometryDataParaview(const Time *, const Geometry *);
static int InitializeTransientParaviewDataFile(ParaviewSet *);
static int WriteSteadyParaviewDataFile(const Time *, ParaviewSet *);
static int WriteParaviewVariableFile(const Geometry *, ParaviewSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ReadGeometryData(const Space *space, const Time *time, const Model *model,
        Geometry *geo)
{
    if (0 == time->restart) { /* if non-restart, read input file */
        NonrestartGeometryReader(geo);
    } else { /* if restart, read the geometry file */
        RestartGeometryReader(time, geo);
    }
    ComputeGeometryParameters(space, model, geo);
    return 0;
}
static int NonrestartGeometryReader(Geometry *geo)
{
    ShowInformation("Reading geometry data ...");
    FILE *filePointer = fopen("artracfd.geo", "r");
    if (NULL == filePointer) {
        FatalError("failed to open geometry file: artracfd.geo...");
    }
    /* read and process file line by line */
    char currentLine[500] = {'\0'}; /* store the current read line */
    char fileName[100] = {'\0'}; /* store current geometry file name */
    int entryCount = 0; /* entry count */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "count begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(geo->totalM));
            if (0 == geo->totalM) {
                ++entryCount;
                break;
            }
            geo->list = AssignStorage(geo->totalM, "Polygon");
            continue;
        }
        if (0 == strncmp(currentLine, "sphere begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(geo->sphereM));
            if (0 == geo->sphereM) {
                continue;
            }
            if (geo->totalM < geo->sphereM) {
                geo->sphereM = geo->totalM;
            }
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%s", fileName);
            ReadSphereFile(fileName, geo);
            if (geo->totalM == geo->sphereM) {
                break;
            }
            continue;
        }
        if (0 == strncmp(currentLine, "STL begin", sizeof currentLine)) {
            ++entryCount;
            int m = geo->sphereM + geo->stlM; /* current geometry pointer */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%s", fileName);
            ReadStlFile(fileName, geo->list + m);
            ReadPolygonStatusData(&filePointer, geo->list + m);
            ++geo->stlM; /* point to the next geometry */
            continue;
        }
    }
    fclose(filePointer); /* close current opened file */
    /* Check missing information section in configuration */
    if (2 + geo->stlM != entryCount) {
        FatalError("missing necessary information section");
    }
    ShowInformation("Session End");
    return 0;
}
static int ReadSphereFile(const char *fileName, Geometry *geo)
{
    FILE *filePointer = fopen(fileName, "r");
    if (NULL == filePointer) {
        FatalError("failed to read sphere geometry file ...");
    }
    for (int m = 0; m < geo->sphereM; ++m) {
        ReadPolygonStatusData(&filePointer, geo->list + m);
        geo->list[m].facetN = 0; /* analytical geometry tag */
        geo->list[m].facet = NULL;
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int ReadPolygonStatusData(FILE **filePointerPointer, Polygon *poly)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    char currentLine[500] = {'\0'}; /* store the current read line */
    /* set format specifier according to the type of Real */
    char format[100] = "%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(format, "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g", sizeof format); /* float type */
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, format,
            &(poly->xc), &(poly->yc), &(poly->zc), &(poly->r),
            &(poly->u), &(poly->v), &(poly->w),
            &(poly->fx), &(poly->fy), &(poly->fz),
            &(poly->rho), &(poly->T), &(poly->cf));
    *filePointerPointer = filePointer; /* updated file pointer */
    return 0;
}
static int WritePolygonStatusData(FILE **filePointerPointer, Polygon *poly)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    fprintf(filePointer, "        <!-- %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g -->\n",
            poly->xc, poly->yc, poly->zc, poly->r,
            poly->u, poly->v, poly->w,
            poly->fx, poly->fy, poly->fz,
            poly->rho, poly->T, poly->cf);
    *filePointerPointer = filePointer; /* updated file pointer */
    return 0;
}
static int RestartGeometryReader(const Time *time, Geometry *geo)
{
    ShowInformation("Restore geometry data ...");
    ReadGeometryDataParaview(time, geo);
    ShowInformation("Session End");
    return 0;
}
static int ReadGeometryDataParaview(const Time *time, Geometry *geo)
{
    ParaviewSet paraSet = { /* initialize ParaviewSet environment */
        .rootName = "geometry", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .intType = "Int32", /* paraview data type */
        .floatType = "Float32", /* paraview data type */
        .byteOrder = "LittleEndian" /* byte order of data */
    };
    snprintf(paraSet.fileName, sizeof(ParaviewString), "%s%05d.vtp", 
            paraSet.rootName, time->restart); 
    FILE *filePointer = fopen(paraSet.fileName, "r");
    if (NULL == filePointer) {
        FatalError("failed to restore geometry file...");
    }
    /* get rid of redundant lines */
    char currentLine[200] = {'\0'}; /* store current line */
    fgets(currentLine, sizeof currentLine, filePointer);
    /* get number of geometries */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %d", &(geo->totalM)); 
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %d", &(geo->sphereM)); 
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %d", &(geo->stlM)); 
    if (0 == geo->totalM) {
        fclose(filePointer);
        return 0;
    }
    geo->list = AssignStorage(geo->totalM, "Polygon");
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 != strncmp(currentLine, "<!-- appended data begin -->", sizeof currentLine)) {
            continue;
        }
        break;
    }
    for (int m = 0; m < geo->sphereM; ++m) {
        ReadPolygonStatusData(&filePointer, geo->list + m);
        geo->list[m].facetN = 0; /* analytical geometry tag */
        geo->list[m].facet = NULL;
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    for (int m = geo->sphereM; m < geo->totalM; ++m) {
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, "%*s %*s NumberOfPolys=\"%d\"", &(geo->list[m].facetN)); 
        geo->list[m].facet = AssignStorage(geo->list[m].facetN, "Facet");
        fgets(currentLine, sizeof currentLine, filePointer);
        fgets(currentLine, sizeof currentLine, filePointer);
        fgets(currentLine, sizeof currentLine, filePointer);
        fgets(currentLine, sizeof currentLine, filePointer);
        fgets(currentLine, sizeof currentLine, filePointer);
        fgets(currentLine, sizeof currentLine, filePointer);
        fgets(currentLine, sizeof currentLine, filePointer);
        ParaviewReal data = 0.0; /* the Paraview data format */
        /* set format specifier according to the type of Real */
        char format[5] = "%lg"; /* default is double type */
        if (sizeof(ParaviewReal) == sizeof(float)) {
            strncpy(format, "%g", sizeof format); /* float type */
        }
        for (int n = 0; n < geo->list[m].facetN; ++m) {
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].x1 = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].y1 = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].z1 = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].x2 = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].y2 = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].z2 = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].x3 = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].y3 = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].z3 = data;
        }
        while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
            CommandLineProcessor(currentLine); /* process current line */
            if (0 != strncmp(currentLine, "<!-- appended data begin -->", sizeof currentLine)) {
                continue;
            }
            break;
        }
        ReadPolygonStatusData(&filePointer, geo->list + m);
        fgets(currentLine, sizeof currentLine, filePointer);
        fgets(currentLine, sizeof currentLine, filePointer);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
int WriteGeometryData(const Time *time, const Geometry *geo)
{
    WriteGeometryDataParaview(time, geo);
    return 0;
}
static int WriteGeometryDataParaview(const Time *time, const Geometry *geo)
{
    ParaviewSet paraSet = { /* initialize ParaviewSet environment */
        .rootName = "geometry", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .intType = "Int32", /* paraview data type */
        .floatType = "Float32", /* paraview data type */
        .byteOrder = "LittleEndian" /* byte order of data */
    };
    snprintf(paraSet.baseName, sizeof(ParaviewString), "%s%05d", 
            paraSet.rootName, time->outputCount); 
    if (0 == time->stepCount) { /* this is the initialization step */
        InitializeTransientParaviewDataFile(&paraSet);
    }
    WriteSteadyParaviewDataFile(time, &paraSet);
    WriteParaviewVariableFile(geo, &paraSet);
    return 0;
}
static int InitializeTransientParaviewDataFile(ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.pvd", 
            paraSet->rootName); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to transient geometry file...");
    }
    /* output information to file */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <Collection>\n");
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteSteadyParaviewDataFile(const Time *time, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.pvd", 
            paraSet->baseName); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to steady geometry file...");
    }
    /* output information to file */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <Collection>\n");
    fprintf(filePointer, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", time->now);
    fprintf(filePointer, "             file=\"%s.vtp\"/>\n", paraSet->baseName);
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "  <!-- Order %d -->\n", time->outputCount);
    fprintf(filePointer, "  <!-- Time %.6g -->\n", time->now);
    fprintf(filePointer, "  <!-- Step %d -->\n", time->stepCount);
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    /*
     * Add the current export to the transient case
     */
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.pvd", 
            paraSet->rootName); 
    filePointer = fopen(paraSet->fileName, "r+");
    if (NULL == filePointer) {
        FatalError("failed to add data to transient geometry file...");
    }
    /* seek the target line for adding information */
    char currentLine[200] = {'\0'}; /* store the current read line */
    int targetLine = 1;
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "</Collection>", sizeof currentLine)) {
            break;
        }
        ++targetLine;
    }
    /* redirect to the target line */
    rewind(filePointer); /* seek to the beginning of the file */
    for (int line = 1; line < targetLine; ++line) {
        fgets(currentLine, sizeof currentLine, filePointer);
    }
    /* append informatiom */
    fprintf(filePointer, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", time->now);
    fprintf(filePointer, "             file=\"%s.vtp\"/>\n", paraSet->baseName);
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteParaviewVariableFile(const Geometry *geo, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.vtp", paraSet->baseName); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to write data file...");
    }
    ParaviewReal data = 0.0; /* paraview scalar data */
    ParaviewReal vector[3] = {0.0}; /* paraview vector data elements */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<!-- M %d -->\n", geo->totalM);
    fprintf(filePointer, "<!-- sphereM %d -->\n", geo->sphereM);
    fprintf(filePointer, "<!-- stlM %d -->\n", geo->stlM);
    fprintf(filePointer, "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <PolyData>");
    fprintf(filePointer, "    <Piece NumberOfPoints=\"%d\" NumberOfPolys=\"%d\">", geo->sphereM, 0);
    fprintf(filePointer, "      <PointData Scalars=\"r\" Vectors=\"Vel\">\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"r\" format=\"ascii\">\n", paraSet->floatType);
    fprintf(filePointer, "          ");
    for (int m = 0; m < geo->sphereM; ++m) {
        data = geo->list[m].r;
        fprintf(filePointer, "%.6g ", data);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"Vel\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(filePointer, "          ");
    for (int m = 0; m < geo->sphereM; ++m) {
        vector[0] = geo->list[m].u;
        vector[1] = geo->list[m].v;
        vector[2] = geo->list[m].w;
        fprintf(filePointer, "%.6g %.6g %.6g ", vector[0], vector[1], vector[2]);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </PointData>\n");
    fprintf(filePointer, "      <CellData>\n");
    fprintf(filePointer, "      </CellData>\n");
    fprintf(filePointer, "      <Points>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"points\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(filePointer, "          ");
    for (int m = 0; m < geo->sphereM; ++m) {
        vector[0] = geo->list[m].xc;
        vector[1] = geo->list[m].yc;
        vector[2] = geo->list[m].zc;
        fprintf(filePointer, "%.6g %.6g %.6g ", vector[0], vector[1], vector[2]);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </Points>\n");
    fprintf(filePointer, "      <Polys>\n");
    fprintf(filePointer, "      </Polys>\n");
    fprintf(filePointer, "      <!-- appended data begin -->\n");
    for (int m = 0; m < geo->sphereM; ++m) {
        WritePolygonStatusData(&filePointer, geo->list + m);
    }
    fprintf(filePointer, "      <!-- appended data end -->\n");
    fprintf(filePointer, "    </Piece>\n");
    for (int m = geo->sphereM; m < geo->totalM; ++m) {
        fprintf(filePointer, "    <Piece NumberOfPoints=\"%d\" NumberOfPolys=\"%d\">", geo->list[m].facetN * 3, geo->list[m].facetN);
        fprintf(filePointer, "      <PointData>\n");
        fprintf(filePointer, "      </PointData>\n");
        fprintf(filePointer, "      <CellData>\n");
        fprintf(filePointer, "      </CellData>\n");
        fprintf(filePointer, "      <Points>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"points\"\n", paraSet->floatType);
        fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
        fprintf(filePointer, "          ");
        for (int n = 0; n < geo->list[m].facetN; ++m) {
            vector[0] = geo->list[m].facet[n].x1;
            vector[1] = geo->list[m].facet[n].y1;
            vector[2] = geo->list[m].facet[n].z1;
            fprintf(filePointer, "%.6g %.6g %.6g ", vector[0], vector[1], vector[2]);
            vector[0] = geo->list[m].facet[n].x2;
            vector[1] = geo->list[m].facet[n].y2;
            vector[2] = geo->list[m].facet[n].z2;
            fprintf(filePointer, "%.6g %.6g %.6g ", vector[0], vector[1], vector[2]);
            vector[0] = geo->list[m].facet[n].x3;
            vector[1] = geo->list[m].facet[n].y3;
            vector[2] = geo->list[m].facet[n].z3;
            fprintf(filePointer, "%.6g %.6g %.6g ", vector[0], vector[1], vector[2]);
        }
        fprintf(filePointer, "\n        </DataArray>\n");
        fprintf(filePointer, "      </Points>\n");
        fprintf(filePointer, "      <Polys>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"ascii\">\n", paraSet->intType);
        fprintf(filePointer, "          ");
        for (int n = 0; n < geo->list[m].facetN; ++m) {
            fprintf(filePointer, "%d %d %d ", 3 * n, 3 * n + 1, 3 * n + 2);
        }
        fprintf(filePointer, "\n        </DataArray>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"ascii\">\n", paraSet->intType);
        fprintf(filePointer, "          ");
        for (int n = 0; n < geo->list[m].facetN; ++m) {
            fprintf(filePointer, "%d ", 3 * (n + 1));
        }
        fprintf(filePointer, "\n        </DataArray>\n");
        fprintf(filePointer, "      </Polys>\n");
        fprintf(filePointer, "      <!-- appended data begin -->\n");
        WritePolygonStatusData(&filePointer, geo->list + m);
        fprintf(filePointer, "      <!-- appended data end -->\n");
        fprintf(filePointer, "    </Piece>\n");
    }
    fprintf(filePointer, "  </PolyData>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

