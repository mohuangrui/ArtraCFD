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
#include "computational_geometry.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void ReadCaseFile(Time *, PvSet *);
static void ReadStructuredData(Space *, const Model *, PvSet *);
static void PointPolyDataReader(const Time *, Geometry *const);
static void ReadPointPolyData(const int, const int, Geometry *const, PvSet *);
static void PolygonPolyDataReader(const Time *, Geometry *const);
static void ReadPolygonPolyData(const int, const int, Geometry *const, PvSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void ReadStructuredDataParaview(Time *time, Space *space, const Model *model)
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
        .scaN = 5,
        .sca = {"rho", "u", "v", "w", "p"},
        .vecN = 0,
        .vec = {{'\0'}},
    };
    snprintf(pvSet.bname, sizeof(PvStr), pvSet.fmt, pvSet.rname, time->dataC);
    ReadCaseFile(time, &pvSet);
    ReadStructuredData(space, model, &pvSet);
    return;
}
static void ReadCaseFile(Time *time, PvSet *pvSet)
{
    snprintf(pvSet->fname, sizeof(PvStr), "%s.pvd", pvSet->bname);
    FILE *fp = Fopen(pvSet->fname, "r");
    ReadInLine(fp, "<!--");
    Sread(fp, 1, ParseFormat("%*s %lg"), &(time->now));
    Sread(fp, 1, "%*s %d", &(time->stepC));
    fclose(fp);
    return;
}
static void ReadStructuredData(Space *space, const Model *model, PvSet *pvSet)
{
    snprintf(pvSet->fname, sizeof(PvStr), "%s%s", pvSet->bname, pvSet->fext);
    FILE *fp = Fopen(pvSet->fname, "r");
    PvReal data = 0.0; /* paraview scalar data */
    const char *fmtI = ParseFormat("%lg");
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    /* get rid of redundant lines */
    ReadInLine(fp, "<PointData>");
    for (int s = 0; s < pvSet->scaN; ++s) {
        Sread(fp, 0, "");
        for (int k = part->ns[PAL][Z][MIN]; k < part->ns[PAL][Z][MAX]; ++k) {
            for (int j = part->ns[PAL][Y][MIN]; j < part->ns[PAL][Y][MAX]; ++j) {
                for (int i = part->ns[PAL][X][MIN]; i < part->ns[PAL][X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    if (0 == s) {
                        /* geometric field initializer */
                        node[idx].did = NONE;
                        node[idx].fid = NONE;
                        node[idx].lid = NONE;
                        node[idx].gst = NONE;
                        memset(node[idx].U, 1, DIMT * sizeof(*node[idx].U));
                        if (InPartBox(k, j, i, part->ns[PIN])) {
                            node[idx].did = 0;
                            node[idx].fid = 0;
                            node[idx].lid = 0;
                            node[idx].gst = 0;
                        }
                    }
                    if (!InPartBox(k, j, i, part->ns[PIO])) {
                        continue;
                    }
                    /* data field initializer */
                    U = node[idx].U[TO];
                    Fscanf(fp, 1, fmtI, &data);
                    switch (s) {
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
        Sread(fp, 0, ""); /* get rid of the end of line of data */
        Sread(fp, 0, "");
    }
    fclose(fp);
    return;
}
void ReadPolyDataParaview(const Time *time, Geometry *const geo)
{
    if (0 != geo->sphN) {
        PointPolyDataReader(time, geo);
    }
    if (0 != geo->stlN) {
        PolygonPolyDataReader(time, geo);
    }
    return;
}
static void PointPolyDataReader(const Time *time, Geometry *const geo)
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
        .scaN = 0,
        .sca = {{'\0'}},
        .vecN = 0,
        .vec = {{'\0'}},
    };
    snprintf(pvSet.bname, sizeof(PvStr), pvSet.fmt, pvSet.rname, time->dataC);
    ReadPointPolyData(0, geo->sphN, geo, &pvSet);
    return;
}
static void ReadPointPolyData(const int pm, const int pn, Geometry *const geo, PvSet *pvSet)
{
    snprintf(pvSet->fname, sizeof(PvStr), "%s%s", pvSet->bname, pvSet->fext);
    FILE *fp = Fopen(pvSet->fname, "r");
    ReadInLine(fp, "<!--");
    ReadPolyStateData(pm, pn, fp, geo);
    fclose(fp);
    return;
}
static void PolygonPolyDataReader(const Time *time, Geometry *const geo)
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
    ReadPolygonPolyData(geo->sphN, geo->totN, geo, &pvSet);
    return;
}
static void ReadPolygonPolyData(const int pm, const int pn, Geometry *const geo, PvSet *pvSet)
{
    snprintf(pvSet->fname, sizeof(PvStr), "%s%s", pvSet->bname, pvSet->fext);
    FILE *fp = Fopen(pvSet->fname, "r");
    PvReal Vec[3] = {0.0}; /* paraview vector data */
    Polyhedron *poly  = NULL;
    const char *fmtJ = ParseFormat("%lg %lg %lg");
    /* get rid of redundant lines */
    ReadInLine(fp, "<PolyData>");
    for (int m = pm; m < pn; ++m) {
        poly = geo->poly + m;
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 1, "%*s %*s %d", &(poly->vertN));
        Sread(fp, 1, "%*s %*s %d", &(poly->edgeN));
        Sread(fp, 1, "%*s %*s %d", &(poly->faceN));
        Sread(fp, 0, "");
        AllocatePolyhedronMemory(poly->vertN, poly->edgeN, poly->faceN, poly);
        poly->edgeN = 0; /* reset edge count before applying edge adding */
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        for (int n = 0; n < poly->vertN; ++n) {
            Fscanf(fp, 3, fmtJ, &(Vec[X]), &(Vec[Y]), &(Vec[Z]));
            poly->v[n][X] = Vec[X];
            poly->v[n][Y] = Vec[Y];
            poly->v[n][Z] = Vec[Z];
        }
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        Sread(fp, 0, "");
        for (int n = 0; n < poly->faceN; ++n) {
            Fscanf(fp, 3, "%d %d %d", &(poly->f[n][0]), &(poly->f[n][1]), &(poly->f[n][2]));
            AddEdge(poly->f[n][0], poly->f[n][1], n, poly);
            AddEdge(poly->f[n][1], poly->f[n][2], n, poly);
            AddEdge(poly->f[n][2], poly->f[n][0], n, poly);
        }
        QuickSortEdge(poly->edgeN, poly->e);
        ReadInLine(fp, "</Piece>");
    }
    ReadInLine(fp, "<!--");
    ReadPolyStateData(pm, pn, fp, geo);
    fclose(fp);
    return;
}
/* a good practice: end file with a newline */

