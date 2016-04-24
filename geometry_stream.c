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
#include "geometry_stream.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "stl.h"
#include "paraview.h"
#include "computational_geometry.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NonrestartGeometryReader(Geometry *);
static int RestartGeometryReader(const Time *, Geometry *);
static int ReadSphereFile(const char *, Geometry *);
static int ReadPolyhedronStatusData(FILE **, Polyhedron *);
static int WritePolyhedronStatusData(FILE **, Polyhedron *);
static int ReadGeometryDataParaview(const Time *, Geometry *);
static int WriteGeometryDataParaview(const Time *, const Geometry *);
static int InitializeTransientParaviewDataFile(ParaviewSet *);
static int WriteSteadyParaviewDataFile(const Time *, ParaviewSet *);
static int WriteParaviewVariableFile(const Geometry *, ParaviewSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ReadGeometryData(Space *space, const Time *time)
{
    if (0 == time->restart) {
        NonrestartGeometryReader(&(space->geo));
    } else {
        RestartGeometryReader(time, &(space->geo));
    }
    ComputeGeometryParameters(space, &(space->geo));
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
    String currentLine = {'\0'}; /* store the current read line */
    String fileName = {'\0'}; /* store current geometry file name */
    int entryCount = 0; /* entry count */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "count begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(geo->totalN));
            if (0 >= geo->totalN) {
                geo->totalN = 0;
                ++entryCount;
                break;
            }
            geo->list = AssignStorage(geo->totalN, "Polyhedron");
            continue;
        }
        if (0 == strncmp(currentLine, "sphere begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(geo->sphereN));
            if (0 >= geo->sphereN) {
                geo->sphereN = 0;
                continue;
            }
            if (geo->totalN < geo->sphereN) {
                geo->sphereN = geo->totalN;
            }
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%s", fileName);
            ReadSphereFile(fileName, geo);
            if (geo->totalN == geo->sphereN) {
                break;
            }
            continue;
        }
        if (0 == strncmp(currentLine, "STL begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%s", fileName);
            ReadStlFile(fileName, geo->list + geo->sphereN + geo->stlN);
            ReadPolyhedronStatusData(&filePointer, geo->list + geo->sphereN + geo->stlN);
            ++geo->stlN; /* point to the next geometry */
            continue;
        }
    }
    fclose(filePointer); /* close current opened file */
    /* Check missing information section in configuration */
    if (2 + geo->stlN != entryCount) {
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
    for (int n = 0; n < geo->sphereN; ++n) {
        ReadPolyhedronStatusData(&filePointer, geo->list + n);
        geo->list[n].facetN = 0; /* analytical geometry tag */
        geo->list[n].facet = NULL;
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int ReadPolyhedronStatusData(FILE **filePointerPointer, Polyhedron *poly)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    String currentLine = {'\0'}; /* store the current read line */
    /* set format specifier according to the type of Real */
    char format[100] = "%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %d";
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(format, "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %d", sizeof format);
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, format,
            &(poly->O[X]), &(poly->O[Y]), &(poly->O[Z]), &(poly->r),
            &(poly->V[X]), &(poly->V[Y]), &(poly->V[Z]),
            &(poly->F[X]), &(poly->F[Y]), &(poly->F[Z]),
            &(poly->rho), &(poly->T), &(poly->cf),
            &(poly->area), &(poly->volume), &(poly->matID));
    *filePointerPointer = filePointer; /* updated file pointer */
    return 0;
}
static int WritePolyhedronStatusData(FILE **filePointerPointer, Polyhedron *poly)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    fprintf(filePointer, "        %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %d\n",
            poly->O[X], poly->O[Y], poly->O[Z], poly->r,
            poly->V[X], poly->V[Y], poly->V[Z],
            poly->F[X], poly->F[Y], poly->F[Z],
            poly->rho, poly->T, poly->cf,
            poly->area, poly->volume, poly->matID);
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
    String currentLine = {'\0'}; /* store current line */
    fgets(currentLine, sizeof currentLine, filePointer);
    /* get number of geometries */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %d", &(geo->totalN)); 
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %d", &(geo->sphereN)); 
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %d", &(geo->stlN)); 
    if (0 == geo->totalN) {
        fclose(filePointer);
        return 0;
    }
    geo->list = AssignStorage(geo->totalN, "Polyhedron");
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 != strncmp(currentLine, "<!-- appended data begin -->", sizeof currentLine)) {
            continue;
        }
        break;
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    for (int m = 0; m < geo->sphereN; ++m) {
        ReadPolyhedronStatusData(&filePointer, geo->list + m);
        geo->list[m].facetN = 0; /* analytical geometry tag */
        geo->list[m].facet = NULL;
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    for (int m = geo->sphereN; m < geo->totalN; ++m) {
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
            geo->list[m].facet[n].P1[X] = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].P1[Y] = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].P1[Z] = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].P2[X] = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].P2[Y] = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].P2[Z] = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].P3[X] = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].P3[Y] = data;
            fscanf(filePointer, format, &data);
            geo->list[m].facet[n].P3[Z] = data;
        }
        while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
            CommandLineProcessor(currentLine); /* process current line */
            if (0 != strncmp(currentLine, "<!-- appended data begin -->", sizeof currentLine)) {
                continue;
            }
            break;
        }
        fgets(currentLine, sizeof currentLine, filePointer);
        ReadPolyhedronStatusData(&filePointer, geo->list + m);
        fgets(currentLine, sizeof currentLine, filePointer);
        fgets(currentLine, sizeof currentLine, filePointer);
        fgets(currentLine, sizeof currentLine, filePointer);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
int WriteGeometryData(const Space *space, const Time *time)
{
    WriteGeometryDataParaview(time, &(space->geo));
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
            paraSet.rootName, time->countOutput); 
    if (0 == time->countStep) { /* this is the initialization step */
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
    fprintf(filePointer, "  <!-- Order %d -->\n", time->countOutput);
    fprintf(filePointer, "  <!-- Time %.6g -->\n", time->now);
    fprintf(filePointer, "  <!-- Step %d -->\n", time->countStep);
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
    String currentLine = {'\0'}; /* store the current read line */
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
    ParaviewVector Vec = {0.0}; /* paraview vector data elements */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<!-- N %d -->\n", geo->totalN);
    fprintf(filePointer, "<!-- sphereN %d -->\n", geo->sphereN);
    fprintf(filePointer, "<!-- stlN %d -->\n", geo->stlN);
    if (0 == geo->totalN) {
        fclose(filePointer);
        return 0;
    }
    fprintf(filePointer, "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <PolyData>");
    fprintf(filePointer, "    <Piece NumberOfPoints=\"%d\" NumberOfPolys=\"%d\">", geo->sphereN, 0);
    fprintf(filePointer, "      <PointData Scalars=\"r\" Vectors=\"Vel\">\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"r\" format=\"ascii\">\n", paraSet->floatType);
    fprintf(filePointer, "          ");
    for (int n = 0; n < geo->sphereN; ++n) {
        data = geo->list[n].r;
        fprintf(filePointer, "%.6g ", data);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"Vel\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"%d\" format=\"ascii\">\n", DIMS);
    fprintf(filePointer, "          ");
    for (int n = 0; n < geo->sphereN; ++n) {
        Vec[X] = geo->list[n].V[X];
        Vec[Y] = geo->list[n].V[Y];
        Vec[Z] = geo->list[n].V[Z];
        fprintf(filePointer, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </PointData>\n");
    fprintf(filePointer, "      <CellData>\n");
    fprintf(filePointer, "      </CellData>\n");
    fprintf(filePointer, "      <Points>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"points\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"%d\" format=\"ascii\">\n", DIMS);
    fprintf(filePointer, "          ");
    for (int n = 0; n < geo->sphereN; ++n) {
        Vec[X] = geo->list[n].O[X];
        Vec[Y] = geo->list[n].O[Y];
        Vec[Z] = geo->list[n].O[Z];
        fprintf(filePointer, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </Points>\n");
    fprintf(filePointer, "      <Polys>\n");
    fprintf(filePointer, "      </Polys>\n");
    fprintf(filePointer, "      <!-- appended data begin -->\n");
    fprintf(filePointer, "      <!-- \n");
    for (int n = 0; n < geo->sphereN; ++n) {
        WritePolyhedronStatusData(&filePointer, geo->list + n);
    }
    fprintf(filePointer, "       -->\n");
    fprintf(filePointer, "      <!-- appended data end -->\n");
    fprintf(filePointer, "    </Piece>\n");
    for (int n = geo->sphereN; n < geo->totalN; ++n) {
        fprintf(filePointer, "    <Piece NumberOfPoints=\"%d\" NumberOfPolys=\"%d\">", geo->list[n].facetN * 3, geo->list[n].facetN);
        fprintf(filePointer, "      <PointData>\n");
        fprintf(filePointer, "      </PointData>\n");
        fprintf(filePointer, "      <CellData>\n");
        fprintf(filePointer, "      </CellData>\n");
        fprintf(filePointer, "      <Points>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"points\"\n", paraSet->floatType);
        fprintf(filePointer, "                   NumberOfComponents=\"%d\" format=\"ascii\">\n", DIMS);
        fprintf(filePointer, "          ");
        for (int m = 0; m < geo->list[n].facetN; ++m) {
            Vec[X] = geo->list[n].facet[m].P1[X];
            Vec[Y] = geo->list[n].facet[m].P1[Y];
            Vec[Z] = geo->list[n].facet[m].P1[Z];
            fprintf(filePointer, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
            Vec[X] = geo->list[n].facet[m].P2[X];
            Vec[Y] = geo->list[n].facet[m].P2[Y];
            Vec[Z] = geo->list[n].facet[m].P2[Z];
            fprintf(filePointer, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
            Vec[X] = geo->list[n].facet[m].P3[X];
            Vec[Y] = geo->list[n].facet[m].P3[Y];
            Vec[Z] = geo->list[n].facet[m].P3[Z];
            fprintf(filePointer, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
        }
        fprintf(filePointer, "\n        </DataArray>\n");
        fprintf(filePointer, "      </Points>\n");
        fprintf(filePointer, "      <Polys>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"ascii\">\n", paraSet->intType);
        fprintf(filePointer, "          ");
        for (int m = 0; m < geo->list[n].facetN; ++m) {
            fprintf(filePointer, "%d %d %d ", 3 * m, 3 * m + 1, 3 * m + 2);
        }
        fprintf(filePointer, "\n        </DataArray>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"ascii\">\n", paraSet->intType);
        fprintf(filePointer, "          ");
        for (int m = 0; m < geo->list[n].facetN; ++m) {
            fprintf(filePointer, "%d ", 3 * (m + 1));
        }
        fprintf(filePointer, "\n        </DataArray>\n");
        fprintf(filePointer, "      </Polys>\n");
        fprintf(filePointer, "      <!-- appended data begin -->\n");
        fprintf(filePointer, "      <!-- \n");
        WritePolyhedronStatusData(&filePointer, geo->list + n);
        fprintf(filePointer, "       -->\n");
        fprintf(filePointer, "      <!-- appended data end -->\n");
        fprintf(filePointer, "    </Piece>\n");
    }
    fprintf(filePointer, "  </PolyData>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

