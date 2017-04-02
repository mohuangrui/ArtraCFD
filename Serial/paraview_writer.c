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
static int PointPolyDataWriter(const Time *, const Geometry *);
static int WritePointPolyData(const int, const int, const Geometry *, ParaviewSet *);
static int PolygonPolyDataWriter(const Time *, const Geometry *);
static int WritePolygonPolyData(const int, const int, const Geometry *, ParaviewSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteStructuredDataParaview(const Time *time, const Space *space, const Model *model)
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
    if (0 == time->stepC) { /* this is the initialization step */
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
    fprintf(filePointer, "  Step %d\n", time->stepC);
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
    WriteToLine(filePointer, "</Collection>");
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
    const char scalar[10][5] = {"rho", "u", "v", "w", "p", "T", "gid", "fid", "lid", "gst"};
    const Partition *restrict part = &(space->part);
    const Node *const node = space->node;
    const Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    IntVec nodeCount = {0}; /* i, j, k node number in each part */
    nodeCount[X] = part->ns[PIN][X][MAX] - part->ns[PIN][X][MIN] - 1; 
    nodeCount[Y] = part->ns[PIN][Y][MAX] - part->ns[PIN][Y][MIN] - 1; 
    nodeCount[Z] = part->ns[PIN][Z][MAX] - part->ns[PIN][Z][MIN] - 1; 
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"StructuredGrid\" version=\"1.0\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 
            0, nodeCount[X], 0, nodeCount[Y], 0, nodeCount[Z]);
    fprintf(filePointer, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", 
            0, nodeCount[X], 0, nodeCount[Y], 0, nodeCount[Z]);
    fprintf(filePointer, "      <PointData>\n");
    for (int count = 0; count < 10; ++count) {
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"%s\" format=\"ascii\">\n", 
                paraSet->floatType, scalar[count]);
        fprintf(filePointer, "          ");
        for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
            for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
                for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
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
                        case 7: /* face flag */
                            data = node[idx].fid;
                            break;
                        case 8: /* layer flag */
                            data = node[idx].lid;
                            break;
                        case 9: /* ghost flag */
                            data = node[idx].gst;
                            break;
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
                U = node[idx].U[TO];
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
int WritePolyDataParaview(const Time *time, const Geometry *geo)
{
    if (0 != geo->sphN) {
        PointPolyDataWriter(time, geo);
    }
    if (0 != geo->stlN) {
        PolygonPolyDataWriter(time, geo);
    }
    return 0;
}
static int PointPolyDataWriter(const Time *time, const Geometry *geo)
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
    if (0 == time->stepC) { /* this is the initialization step */
        InitializeTransientCaseFile(&paraSet);
    }
    WriteCaseFile(time, &paraSet);
    WritePointPolyData(0, geo->sphN, geo, &paraSet);
    return 0;
}
static int WritePointPolyData(const int start, const int end, const Geometry *geo, ParaviewSet *paraSet)
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
    fprintf(filePointer, "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"1\" NumberOfPolys=\"0\">\n", (end - start));
    fprintf(filePointer, "      <PointData>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"r\" format=\"ascii\">\n", paraSet->floatType);
    fprintf(filePointer, "          ");
    for (int n = start; n < end; ++n) {
        data = geo->poly[n].r;
        fprintf(filePointer, "%.6g ", data);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"gid\" format=\"ascii\">\n", paraSet->intType);
    fprintf(filePointer, "          ");
    for (int n = start; n < end; ++n) {
        fprintf(filePointer, "%d ", n + 1);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"Vel\"\n", paraSet->floatType);
    fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
    fprintf(filePointer, "          ");
    for (int n = start; n < end; ++n) {
        Vec[X] = geo->poly[n].V[TO][X];
        Vec[Y] = geo->poly[n].V[TO][Y];
        Vec[Z] = geo->poly[n].V[TO][Z];
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
    for (int n = start; n < end; ++n) {
        Vec[X] = geo->poly[n].O[X];
        Vec[Y] = geo->poly[n].O[Y];
        Vec[Z] = geo->poly[n].O[Z];
        fprintf(filePointer, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </Points>\n");
    fprintf(filePointer, "      <Verts>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"ascii\">\n", paraSet->intType);
    fprintf(filePointer, "          ");
    for (int n = start; n < end; ++n) {
        fprintf(filePointer, "%d ", n);
    }
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"ascii\">\n", paraSet->intType);
    fprintf(filePointer, "          ");
    fprintf(filePointer, "%d ", (end - start));
    fprintf(filePointer, "\n        </DataArray>\n");
    fprintf(filePointer, "      </Verts>\n");
    fprintf(filePointer, "      <Polys>\n");
    fprintf(filePointer, "      </Polys>\n");
    fprintf(filePointer, "    </Piece>\n");
    fprintf(filePointer, "  </PolyData>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fprintf(filePointer, "<!--\n");
    WritePolyhedronStateData(start, end, filePointer, geo);
    fprintf(filePointer, "-->\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int PolygonPolyDataWriter(const Time *time, const Geometry *geo)
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
    if (0 == time->stepC) { /* this is the initialization step */
        InitializeTransientCaseFile(&paraSet);
    }
    WriteCaseFile(time, &paraSet);
    WritePolygonPolyData(geo->sphN, geo->totN, geo, &paraSet);
    return 0;
}
static int WritePolygonPolyData(const int start, const int end, const Geometry *geo, ParaviewSet *paraSet)
{
    snprintf(paraSet->fileName, sizeof(ParaviewString), "%s%s", paraSet->baseName, paraSet->fileExt); 
    FILE *filePointer = fopen(paraSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    ParaviewReal Vec[3] = {0.0}; /* paraview vector data */
    const Polyhedron *poly = NULL;
    fprintf(filePointer, "<?xml version=\"1.0\"?>\n");
    fprintf(filePointer, "<VTKFile type=\"PolyData\" version=\"1.0\"\n");
    fprintf(filePointer, "         byte_order=\"%s\">\n", paraSet->byteOrder);
    fprintf(filePointer, "  <PolyData>\n");
    for (int m = start; m < end; ++m) {
        poly = geo->poly + m;
        fprintf(filePointer, "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfPolys=\"%d\">\n",
                poly->vertN, poly->faceN);
        fprintf(filePointer, "      <!--\n");
        fprintf(filePointer, "        vertN = %d\n", poly->vertN);
        fprintf(filePointer, "        edgeN = %d\n", poly->edgeN);
        fprintf(filePointer, "        faceN = %d\n", poly->faceN);
        fprintf(filePointer, "      -->\n");
        fprintf(filePointer, "      <PointData>\n");
        fprintf(filePointer, "      </PointData>\n");
        fprintf(filePointer, "      <CellData>\n");
        fprintf(filePointer, "      </CellData>\n");
        fprintf(filePointer, "      <Points>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"points\"\n", paraSet->floatType);
        fprintf(filePointer, "                   NumberOfComponents=\"3\" format=\"ascii\">\n");
        fprintf(filePointer, "          ");
        for (int n = 0; n < poly->vertN; ++n) {
            Vec[X] = poly->v[n][X];
            Vec[Y] = poly->v[n][Y];
            Vec[Z] = poly->v[n][Z];
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
        for (int n = 0; n < poly->faceN; ++n) {
            fprintf(filePointer, "%d %d %d ", poly->f[n][0], poly->f[n][1], poly->f[n][2]);
        }
        fprintf(filePointer, "\n        </DataArray>\n");
        fprintf(filePointer, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"ascii\">\n", paraSet->intType);
        fprintf(filePointer, "          ");
        for (int n = 0; n < poly->faceN; ++n) {
            fprintf(filePointer, "%d ", 3 * (n + 1));
        }
        fprintf(filePointer, "\n        </DataArray>\n");
        fprintf(filePointer, "      </Polys>\n");
        fprintf(filePointer, "    </Piece>\n");
    }
    fprintf(filePointer, "  </PolyData>\n");
    fprintf(filePointer, "</VTKFile>\n");
    fprintf(filePointer, "<!--\n");
    WritePolyhedronStateData(start, end, filePointer, geo);
    fprintf(filePointer, "-->\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
int WritePolyhedronStateData(const int start, const int end, FILE *filePointer, const Geometry *geo)
{
    const char formatI[100] = "  %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %d\n";
    const char formatJ[100] = "  %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n";
    const Polyhedron *poly = NULL;
    for (int n = start; n < end; ++n) {
        poly = geo->poly + n;
        fprintf(filePointer, formatI,
                poly->O[X], poly->O[Y], poly->O[Z], poly->r,
                poly->V[TO][X], poly->V[TO][Y], poly->V[TO][Z],
                poly->W[TO][X], poly->W[TO][Y], poly->W[TO][Z],
                poly->rho, poly->T, poly->cf,
                poly->area, poly->volume, poly->mid);
        fprintf(filePointer, formatJ,
                poly->at[TO][X], poly->at[TO][Y], poly->at[TO][Z],
                poly->ar[TO][X], poly->ar[TO][Y], poly->ar[TO][Z],
                poly->at[TN][X], poly->at[TN][Y], poly->at[TN][Z],
                poly->g[X], poly->g[Y], poly->g[Z],
                poly->ar[TN][X], poly->ar[TN][Y], poly->ar[TN][Z],
                poly->to);
    }
    return 0;
}
/* a good practice: end file with a newline */

