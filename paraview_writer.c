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
static int InitializeTransientCaseFile(ParaviewSet *);
static int WriteCaseFile(const Time *, ParaviewSet *);
static int WriteStructuredData(const Space *, const Model *, ParaviewSet *);
static int PointPolyDataWriter(const Geometry *, const Time *);
static int WritePointPolyData(const Geometry *, ParaviewSet *);
static int PolygonPolyDataWriter(const Geometry *, const Time *);
static int WritePolygonPolyData(const Geometry *, ParaviewSet *);
static int WriteToLine(FILE **, const char *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteStructuredDataParaview(const Space *space, const Time *time, const Model *model)
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
    if (0 == time->countStep) { /* this is the initialization step */
        InitializeTransientCaseFile(&paraSet);
    }
    WriteCaseFile(time, &paraSet);
    WriteStructuredData(space, model, &paraSet);
    return 0;
}
static int InitializeTransientCaseFile(ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.pvd", 
            paraSet->rootName); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to initialize transient case file...");
    }
    /* output information to file */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"Collection\" version=\"1.0\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <Collection>\n");
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteCaseFile(const Time *time, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.pvd", 
            paraSet->baseName); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to open case file...");
    }
    /* output information to file */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"Collection\" version=\"1.0\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <Collection>\n");
    fprintf(filePointer, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", time->now);
    fprintf(filePointer, "             file=\"%s%s\"/>\n", paraSet->baseName, paraSet->fileExt);
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fprintf(filePointer, "<!--\n");
    fprintf(filePointer, "  Time %.6g\n", time->now);
    fprintf(filePointer, "  Step %d\n", time->countStep);
    fprintf(filePointer, "-->\n");
    fclose(filePointer); /* close current opened file */
    /*
     * Add the current export to the transient case
     */
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s.pvd", 
            paraSet->rootName); 
    filePointer = fopen(paraSet->fileName, "r+");
    if (NULL == filePointer) {
        FatalError("failed to add data to transient case file...");
    }
    /* seek the target line for adding information */
    WriteToLine(&filePointer, "</Collection>");
    /* append informatiom */
    fprintf(filePointer, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", time->now);
    fprintf(filePointer, "             file=\"%s%s\"/>\n", paraSet->baseName, paraSet->fileExt);
    fprintf(filePointer, "  </Collection>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WriteStructuredData(const Space *space, const Model *model, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%s", paraSet->baseName, paraSet->fileExt); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    ParaviewReal data = 0.0; /* paraview scalar data */
    ParaviewReal Vec[3] = {0.0}; /* paraview vector data */
    /* the scalar variables */
    const char scalar[7][5] = {"rho", "u", "v", "w", "p", "T", "geo"};
    int idx = 0; /* linear array index math variable */
    const Node *node = space->node;
    const Real *restrict U = NULL;
    const Partition *restrict part = &(space->part);
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"StructuredGrid\" version=\"1.0\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 
            0, part->m[X] - 2, 0, part->m[Y] - 2, 0, part->m[Z] - 2);
    fprintf(filePointer, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", 
            0, part->m[X] - 2, 0, part->m[Y] - 2, 0, part->m[Z] - 2);
    fprintf(filePointer, "      <PointData Scalars=\"rho\" Vectors=\"Vel\">\n");
    for (int count = 0; count < 7; ++count) {
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"%s\" format=\"ascii\">\n", 
                paraSet->floatType, scalar[count]);
        fprintf(filePointer, "          ");
        for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
            for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
                for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    U = node[idx].U[C];
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
                            data = node[idx].geoID;
                        default:
                            break;
                    }
                    fprintf(filePointer, "%.6g ", data);
                }
            }
        }
        fprintf(filePointer, "\n        </DataArray>\n");
    }
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"Vel\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(filePointer, "          ");
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                U = node[idx].U[C];
                Vec[X] = U[1] / U[0];
                Vec[Y] = U[2] / U[0];
                Vec[Z] = U[3] / U[0];
                fprintf(filePointer, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
            }
        }
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </PointData>\n");
    fprintf(filePointer, "      <CellData>\n");
    fprintf(filePointer, "      </CellData>\n");
    fprintf(filePointer, "      <Points>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"points\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(filePointer, "          ");
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                Vec[X] = PointSpace(i, part->domain[X][MIN], part->d[X], part->ng);
                Vec[Y] = PointSpace(j, part->domain[Y][MIN], part->d[Y], part->ng);
                Vec[Z] = PointSpace(k, part->domain[Z][MIN], part->d[Z], part->ng);
                fprintf(filePointer, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
            }
        }
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </Points>\n");
    fprintf(filePointer, "    </Piece>\n");
    fprintf(filePointer, "  </StructuredGrid>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
int WritePolyDataParaview(const Geometry *geo, const Time *time)
{
    if (0 != geo->sphereN) {
        PointPolyDataWriter(geo, time);
    }
    if (0 != geo->stlN) {
        PolygonPolyDataWriter(geo, time);
    }
    return 0;
}
static int PointPolyDataWriter(const Geometry *geo, const Time *time)
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
    if (0 == time->countStep) { /* this is the initialization step */
        InitializeTransientCaseFile(&paraSet);
    }
    WriteCaseFile(time, &paraSet);
    WritePointPolyData(geo, &paraSet);
    return 0;
}
static int WritePointPolyData(const Geometry *geo, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%s", paraSet->baseName, paraSet->fileExt); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    ParaviewReal data = 0.0; /* paraview scalar data */
    ParaviewReal Vec[3] = {0.0}; /* paraview vector data */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"PolyData\" version=\"1.0\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <PolyData>\n");
    fprintf(filePointer, "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"1\" NumberOfPolys=\"0\">\n", geo->sphereN);
    fprintf(filePointer, "      <PointData Scalars=\"r\" Vectors=\"Vel\">\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"r\" format=\"ascii\">\n", paraSet->floatType);
    fprintf(filePointer, "          ");
    for (int n = 0; n < geo->sphereN; ++n) {
        data = geo->list[n].r;
        fprintf(filePointer, "%.6g ", data);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"Vel\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
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
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(filePointer, "          ");
    for (int n = 0; n < geo->sphereN; ++n) {
        Vec[X] = geo->list[n].O[X];
        Vec[Y] = geo->list[n].O[Y];
        Vec[Z] = geo->list[n].O[Z];
        fprintf(filePointer, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </Points>\n");
    fprintf(filePointer, "      <Verts>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"ascii\">\n", paraSet->intType);
    fprintf(filePointer, "          ");
    for (int n = 0; n < geo->sphereN; ++n) {
        fprintf(filePointer, "%d ", n);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"ascii\">\n", paraSet->intType);
    fprintf(filePointer, "          ");
    fprintf(filePointer, "%d ", geo->sphereN);
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </Verts>\n");
    fprintf(filePointer, "      <Polys>\n");
    fprintf(filePointer, "      </Polys>\n");
    fprintf(filePointer, "    </Piece>\n");
    fprintf(filePointer, "  </PolyData>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fprintf(filePointer, "<!--\n");
    WritePolyhedronStateData(&filePointer, geo, 0, geo->sphereN);
    fprintf(filePointer, "-->\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int PolygonPolyDataWriter(const Geometry *geo, const Time *time)
{
    ParaviewSet paraSet = { /* initialize ParaviewSet environment */
        .rootName = "geo_stl", /* data file root name */
        .baseName = {'\0'}, /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .fileExt = ".vtp", /* data file extension */
        .intType = "Int32", /* paraview int type */
        .floatType = "Float32", /* paraview float type */
        .byteOrder = "LittleEndian" /* byte order of data */
    };
    snprintf(paraSet.baseName, sizeof(ParaviewString), "%s%05d", 
            paraSet.rootName, time->countOutput); 
    if (0 == time->countStep) { /* this is the initialization step */
        InitializeTransientCaseFile(&paraSet);
    }
    WriteCaseFile(time, &paraSet);
    WritePolygonPolyData(geo, &paraSet);
    return 0;
}
static int WritePolygonPolyData(const Geometry *geo, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%s", paraSet->baseName, paraSet->fileExt); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    ParaviewReal data = 0.0; /* paraview scalar data */
    ParaviewReal Vec[3] = {0.0}; /* paraview vector data */
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"PolyData\" version=\"1.0\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <PolyData>\n");
    for (int n = geo->sphereN; n < geo->totalN; ++n) {
        fprintf(filePointer, "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfPolys=\"%d\">\n",
                geo->list[n].facetN * 3, geo->list[n].facetN);
        fprintf(filePointer, "      <PointData>\n");
        fprintf(filePointer, "      </PointData>\n");
        fprintf(filePointer, "      <CellData>\n");
        fprintf(filePointer, "      </CellData>\n");
        fprintf(filePointer, "      <Points>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"points\"\n", paraSet->floatType);
        fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
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
        fprintf(filePointer, "      <Verts>\n");
        fprintf(filePointer, "      </Verts>\n");
        fprintf(filePointer, "      <Polys>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"ascii\">\n",
                paraSet->intType);
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
        fprintf(filePointer, "    </Piece>\n");
    }
    fprintf(filePointer, "  </PolyData>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fprintf(filePointer, "<!--\n");
    WritePolyhedronStateData(&filePointer, geo, geo->sphereN, geo->totalN);
    fprintf(filePointer, "-->\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
int WritePolyhedronStateData(FILE **filePointerPointer, const Geometry *geo, const int start, const int end)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    for (int n = start; n < end; ++n) {
        fprintf(filePointer, "        %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %d\n",
                geo->list[n].O[X], geo->list[n].O[Y], geo->list[n].O[Z], geo->list[n].r,
                geo->list[n].V[X], geo->list[n].V[Y], geo->list[n].V[Z],
                geo->list[n].F[X], geo->list[n].F[Y], geo->list[n].F[Z],
                geo->list[n].rho, geo->list[n].T, geo->list[n].cf,
                geo->list[n].area, geo->list[n].volume, geo->list[n].matID);
    }
    *filePointerPointer = filePointer; /* updated file pointer */
    return 0;
}
static int WriteToLine(FILE **filePointerPointer, const char *lineString)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    String currentLine = {'\0'}; /* store the current read line */
    int targetLine = 1;
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, lineString, sizeof currentLine)) {
            break;
        }
        ++targetLine;
    }
    /* redirect to the target line */
    filePointer = *filePointerPointer; /* reset */
    for (int count = 1; count < targetLine; ++count) {
        fgets(currentLine, sizeof currentLine, filePointer);
    }
    *filePointerPointer = filePointer; /* updated file pointer */
    return 0;
}
/* a good practice: end file with a newline */

