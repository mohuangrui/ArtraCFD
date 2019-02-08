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
#include "data_stream.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void InitializeTransientCaseFile(PvSet *);
static void WriteCaseFile(const Time *, PvSet *);
static void WriteStructuredData(const Space *, const Model *, PvSet *);
static void PointPolyDataWriter(const Time *, const Geometry *const);
static void WritePointPolyData(const int, const int, const Geometry *const, PvSet *);
static void PolygonPolyDataWriter(const Time *, const Geometry *const);
static void WritePolygonPolyData(const int, const int, const Geometry *const, PvSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void WriteStructuredDataParaview(const Time *time, const Space *space, const Model *model)
{
    PvSet pvSet = { /* initialize environment */
        .rname = "field",
        .bname = {'\0'},
        .fname = {'\0'},
        .fext = ".vts",
        .fmt = "%s%05d",
        .intType = "Int32",
        .floatType = "Float32",
        .byteOrder = "LittleEndian",
        .scaN = 10,
        .sca = {"rho", "u", "v", "w", "p", "T", "did", "fid", "lid", "gst"},
        .vecN = 1,
        .vec = {"Vel"},
    };
    snprintf(pvSet.bname, sizeof(PvStr), pvSet.fmt, pvSet.rname, time->dataC);
    if (0 == time->stepC) { /* initialization step */
        InitializeTransientCaseFile(&pvSet);
    }
    WriteCaseFile(time, &pvSet);
    WriteStructuredData(space, model, &pvSet);
    return;
}
static void InitializeTransientCaseFile(PvSet *pvSet)
{
    snprintf(pvSet->fname, sizeof(PvStr), "%s.pvd", pvSet->rname);
    FILE *fp = Fopen(pvSet->fname, "w");
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"%s\">\n", pvSet->byteOrder);
    fprintf(fp, "  <Collection>\n");
    fprintf(fp, "  </Collection>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
    return;
}
static void WriteCaseFile(const Time *time, PvSet *pvSet)
{
    snprintf(pvSet->fname, sizeof(PvStr), "%s.pvd", pvSet->bname);
    FILE *fp = Fopen(pvSet->fname, "w");
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"%s\">\n", pvSet->byteOrder);
    fprintf(fp, "  <Collection>\n");
    fprintf(fp, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", time->now);
    fprintf(fp, "             file=\"%s%s\"/>\n", pvSet->bname, pvSet->fext);
    fprintf(fp, "  </Collection>\n");
    fprintf(fp, "</VTKFile>\n");
    fprintf(fp, "<!--\n");
    fprintf(fp, "  Time %.6g\n", time->now);
    fprintf(fp, "  Step %d\n", time->stepC);
    fprintf(fp, "-->\n");
    fclose(fp);
    /* add case to the transient case */
    snprintf(pvSet->fname, sizeof(PvStr), "%s.pvd", pvSet->rname);
    fp = Fopen(pvSet->fname, "r+");
    /* seek the target line for adding information */
    WriteToLine(fp, "</Collection>");
    /* append informatiom */
    fprintf(fp, "    <DataSet timestep=\"%.6g\" group=\"\" part=\"0\"\n", time->now);
    fprintf(fp, "             file=\"%s%s\"/>\n", pvSet->bname, pvSet->fext);
    fprintf(fp, "  </Collection>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
    return;
}
static void WriteStructuredData(const Space *space, const Model *model, PvSet *pvSet)
{
    snprintf(pvSet->fname, sizeof(PvStr), "%s%s", pvSet->bname, pvSet->fext);
    FILE *fp = Fopen(pvSet->fname, "w");
    PvReal data = 0.0; /* paraview scalar data */
    PvReal Vec[3] = {0.0}; /* paraview vector data */
    const Partition *const part = &(space->part);
    const Node *const node = space->node;
    const Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    IntVec ne = {0}; /* i, j, k node number in each part */
    ne[X] = part->ns[PIO][X][MAX] - part->ns[PIO][X][MIN] - 1;
    ne[Y] = part->ns[PIO][Y][MAX] - part->ns[PIO][Y][MIN] - 1;
    ne[Z] = part->ns[PIO][Z][MAX] - part->ns[PIO][Z][MIN] - 1;
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"%s\">\n", pvSet->byteOrder);
    fprintf(fp, "  <StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n", 0, ne[X], 0, ne[Y], 0, ne[Z]);
    fprintf(fp, "    <Piece Extent=\"%d %d %d %d %d %d\">\n", 0, ne[X], 0, ne[Y], 0, ne[Z]);
    fprintf(fp, "      <PointData>\n");
    for (int s = 0; s < pvSet->scaN; ++s) {
        fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" format=\"ascii\">\n", pvSet->floatType, pvSet->sca[s]);
        fprintf(fp, "          ");
        for (int k = part->ns[PIO][Z][MIN]; k < part->ns[PIO][Z][MAX]; ++k) {
            for (int j = part->ns[PIO][Y][MIN]; j < part->ns[PIO][Y][MAX]; ++j) {
                for (int i = part->ns[PIO][X][MIN]; i < part->ns[PIO][X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    U = node[idx].U[TO];
                    switch (s) {
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
                            data = node[idx].did;
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
                    fprintf(fp, "%.6g ", data);
                }
            }
        }
        fprintf(fp, "\n        </DataArray>\n");
    }
    for (int s = 0; s < pvSet->vecN; ++s) {
        fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", pvSet->floatType, pvSet->vec[s]);
        fprintf(fp, "          ");
        for (int k = part->ns[PIO][Z][MIN]; k < part->ns[PIO][Z][MAX]; ++k) {
            for (int j = part->ns[PIO][Y][MIN]; j < part->ns[PIO][Y][MAX]; ++j) {
                for (int i = part->ns[PIO][X][MIN]; i < part->ns[PIO][X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    U = node[idx].U[TO];
                    Vec[X] = U[1] / U[0];
                    Vec[Y] = U[2] / U[0];
                    Vec[Z] = U[3] / U[0];
                    fprintf(fp, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
                }
            }
        }
        fprintf(fp, "\n        </DataArray>\n");
    }
    fprintf(fp, "      </PointData>\n");
    fprintf(fp, "      <CellData>\n");
    fprintf(fp, "      </CellData>\n");
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\">\n", pvSet->floatType);
    fprintf(fp, "          ");
    for (int k = part->ns[PIO][Z][MIN]; k < part->ns[PIO][Z][MAX]; ++k) {
        for (int j = part->ns[PIO][Y][MIN]; j < part->ns[PIO][Y][MAX]; ++j) {
            for (int i = part->ns[PIO][X][MIN]; i < part->ns[PIO][X][MAX]; ++i) {
                Vec[X] = MapPoint(i, part->domain[X][MIN], part->d[X], part->ng[X]);
                Vec[Y] = MapPoint(j, part->domain[Y][MIN], part->d[Y], part->ng[Y]);
                Vec[Z] = MapPoint(k, part->domain[Z][MIN], part->d[Z], part->ng[Z]);
                fprintf(fp, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
            }
        }
    }
    fprintf(fp, "\n        </DataArray>\n");
    fprintf(fp, "      </Points>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </StructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
    fclose(fp);
    return;
}
void WritePolyDataParaview(const Time *time, const Geometry *const geo)
{
    if (0 != geo->sphN) {
        PointPolyDataWriter(time, geo);
    }
    if (0 != geo->stlN) {
        PolygonPolyDataWriter(time, geo);
    }
    return;
}
static void PointPolyDataWriter(const Time *time, const Geometry *const geo)
{
    PvSet pvSet = { /* initialize environment */
        .rname = "geo_sph",
        .bname = {'\0'},
        .fname = {'\0'},
        .fext = ".vtp",
        .fmt = "%s%05d",
        .intType = "Int32",
        .floatType = "Float32",
        .byteOrder = "LittleEndian",
        .scaN = 2,
        .sca = {"r", "did"},
        .vecN = 1,
        .vec = {"Vel"},
    };
    snprintf(pvSet.bname, sizeof(PvStr), pvSet.fmt, pvSet.rname, time->dataC);
    if (0 == time->stepC) { /* initialization step */
        InitializeTransientCaseFile(&pvSet);
    }
    WriteCaseFile(time, &pvSet);
    WritePointPolyData(0, geo->sphN, geo, &pvSet);
    return;
}
static void WritePointPolyData(const int pm, const int pn, const Geometry *const geo, PvSet *pvSet)
{
    snprintf(pvSet->fname, sizeof(PvStr), "%s%s", pvSet->bname, pvSet->fext);
    FILE *fp = Fopen(pvSet->fname, "w");
    PvReal data = 0.0; /* paraview scalar data */
    PvReal Vec[3] = {0.0}; /* paraview vector data */
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"%s\">\n", pvSet->byteOrder);
    fprintf(fp, "  <PolyData>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"1\" NumberOfPolys=\"0\">\n", (pn - pm));
    fprintf(fp, "      <PointData>\n");
    for (int s = 0; s < pvSet->scaN; ++s) {
        fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" format=\"ascii\">\n", pvSet->floatType, pvSet->sca[s]);
        fprintf(fp, "          ");
        for (int n = pm; n < pn; ++n) {
            switch (s) {
                case 0:
                    data = geo->poly[n].r;
                    break;
                case 1:
                    data = n + 1;
                    break;
                default:
                    break;
            }
            fprintf(fp, "%.6g ", data);
        }
        fprintf(fp, "\n        </DataArray>\n");
    }
    for (int s = 0; s < pvSet->vecN; ++s) {
        fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"3\" format=\"ascii\">\n", pvSet->floatType, pvSet->vec[s]);
        fprintf(fp, "          ");
        for (int n = pm; n < pn; ++n) {
            Vec[X] = geo->poly[n].V[TO][X];
            Vec[Y] = geo->poly[n].V[TO][Y];
            Vec[Z] = geo->poly[n].V[TO][Z];
            fprintf(fp, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
        }
        fprintf(fp, "\n        </DataArray>\n");
    }
    fprintf(fp, "      </PointData>\n");
    fprintf(fp, "      <CellData>\n");
    fprintf(fp, "      </CellData>\n");
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\">\n", pvSet->floatType);
    fprintf(fp, "          ");
    for (int n = pm; n < pn; ++n) {
        Vec[X] = geo->poly[n].O[X];
        Vec[Y] = geo->poly[n].O[Y];
        Vec[Z] = geo->poly[n].O[Z];
        fprintf(fp, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
    }
    fprintf(fp, "\n        </DataArray>\n");
    fprintf(fp, "      </Points>\n");
    fprintf(fp, "      <Verts>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"ascii\">\n", pvSet->intType);
    fprintf(fp, "          ");
    for (int n = pm; n < pn; ++n) {
        fprintf(fp, "%d ", n);
    }
    fprintf(fp, "\n        </DataArray>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"ascii\">\n", pvSet->intType);
    fprintf(fp, "          ");
    fprintf(fp, "%d ", (pn - pm));
    fprintf(fp, "\n        </DataArray>\n");
    fprintf(fp, "      </Verts>\n");
    fprintf(fp, "      <Polys>\n");
    fprintf(fp, "      </Polys>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </PolyData>\n");
    fprintf(fp, "</VTKFile>\n");
    fprintf(fp, "<!--\n");
    WritePolyStateData(pm, pn, fp, geo);
    fprintf(fp, "-->\n");
    fclose(fp);
    return;
}
static void PolygonPolyDataWriter(const Time *time, const Geometry *const geo)
{
    PvSet pvSet = { /* initialize environment */
        .rname = "geo_stl",
        .bname = {'\0'},
        .fname = {'\0'},
        .fext = ".vtp",
        .fmt = "%s%05d",
        .intType = "Int32",
        .floatType = "Float32",
        .byteOrder = "LittleEndian",
        .scaN = 0,
        .sca = {{'\0'}},
        .vecN = 0,
        .vec = {{'\0'}},
    };
    snprintf(pvSet.bname, sizeof(PvStr), pvSet.fmt, pvSet.rname, time->dataC);
    if (0 == time->stepC) { /* initialization step */
        InitializeTransientCaseFile(&pvSet);
    }
    WriteCaseFile(time, &pvSet);
    WritePolygonPolyData(geo->sphN, geo->totN, geo, &pvSet);
    return;
}
static void WritePolygonPolyData(const int pm, const int pn, const Geometry *const geo, PvSet *pvSet)
{
    snprintf(pvSet->fname, sizeof(PvStr), "%s%s", pvSet->bname, pvSet->fext);
    FILE *fp = Fopen(pvSet->fname, "w");
    PvReal Vec[3] = {0.0}; /* paraview vector data */
    const Polyhedron *poly = NULL;
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"%s\">\n", pvSet->byteOrder);
    fprintf(fp, "  <PolyData>\n");
    for (int m = pm; m < pn; ++m) {
        poly = geo->poly + m;
        fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" NumberOfPolys=\"%d\">\n", poly->vertN, poly->faceN);
        fprintf(fp, "      <!--\n");
        fprintf(fp, "        vertN = %d\n", poly->vertN);
        fprintf(fp, "        edgeN = %d\n", poly->edgeN);
        fprintf(fp, "        faceN = %d\n", poly->faceN);
        fprintf(fp, "      -->\n");
        fprintf(fp, "      <PointData>\n");
        fprintf(fp, "      </PointData>\n");
        fprintf(fp, "      <CellData>\n");
        fprintf(fp, "      </CellData>\n");
        fprintf(fp, "      <Points>\n");
        fprintf(fp, "        <DataArray type=\"%s\" Name=\"points\" NumberOfComponents=\"3\" format=\"ascii\">\n", pvSet->floatType);
        fprintf(fp, "          ");
        for (int n = 0; n < poly->vertN; ++n) {
            Vec[X] = poly->v[n][X];
            Vec[Y] = poly->v[n][Y];
            Vec[Z] = poly->v[n][Z];
            fprintf(fp, "%.6g %.6g %.6g ", Vec[X], Vec[Y], Vec[Z]);
        }
        fprintf(fp, "\n        </DataArray>\n");
        fprintf(fp, "      </Points>\n");
        fprintf(fp, "      <Verts>\n");
        fprintf(fp, "      </Verts>\n");
        fprintf(fp, "      <Polys>\n");
        fprintf(fp, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"ascii\">\n", pvSet->intType);
        fprintf(fp, "          ");
        for (int n = 0; n < poly->faceN; ++n) {
            fprintf(fp, "%d %d %d ", poly->f[n][0], poly->f[n][1], poly->f[n][2]);
        }
        fprintf(fp, "\n        </DataArray>\n");
        fprintf(fp, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"ascii\">\n", pvSet->intType);
        fprintf(fp, "          ");
        for (int n = 0; n < poly->faceN; ++n) {
            fprintf(fp, "%d ", 3 * (n + 1));
        }
        fprintf(fp, "\n        </DataArray>\n");
        fprintf(fp, "      </Polys>\n");
        fprintf(fp, "    </Piece>\n");
    }
    fprintf(fp, "  </PolyData>\n");
    fprintf(fp, "</VTKFile>\n");
    fprintf(fp, "<!--\n");
    WritePolyStateData(pm, pn, fp, geo);
    fprintf(fp, "-->\n");
    fclose(fp);
    return;
}
/* a good practice: end file with a newline */

