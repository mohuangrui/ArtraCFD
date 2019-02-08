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
#include "data_stream.h"
#include "computational_geometry.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void ReadCaseFile(Time *, EnSet *);
static void ReadStructuredData(Space *, const Model *, EnSet *);
static void PointPolyDataReader(const Time *, Geometry *const);
static void PolygonPolyDataReader(const Time *, Geometry *const);
static void ReadPolygonPolyData(const int, const int, Geometry *const, EnSet *);
static void ReadPolyState(const int, const int, Geometry *const, EnSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void ReadStructuredDataEnsight(Time *time, Space *space, const Model *model)
{
    EnSet enSet = { /* initialize environment */
        .rname = "field",
        .bname = {'\0'},
        .fname = {'\0'},
        .str = {'\0'},
        .fmt = "%s%05d",
        .gtag = {'\0'},
        .vtag = "*****",
        .dtype = "block",
        .part = {PIO, PIO + 1},
        .scaN = 5,
        .sca = {"rho", "u", "v", "w", "p"},
        .vecN = 0,
        .vec = {{'\0'}},
    };
    snprintf(enSet.bname, sizeof(EnStr), enSet.fmt, enSet.rname, time->dataC);
    ReadCaseFile(time, &enSet);
    ReadStructuredData(space, model, &enSet);
    return;
}
static void ReadCaseFile(Time *time, EnSet *enSet)
{
    snprintf(enSet->fname, sizeof(EnStr), "%s.case", enSet->bname);
    FILE *fp = Fopen(enSet->fname, "r");
    ReadInLine(fp, "VARIABLE");
    Sread(fp, 1, ParseFormat("%*s %*s %*s %*s %lg"), &(time->now));
    Sread(fp, 1, "%*s %*s %*s %*s %d", &(time->stepC));
    fclose(fp);
    return;
}
static void ReadStructuredData(Space *space, const Model *model, EnSet *enSet)
{
    FILE *fp = NULL;
    EnReal data = 0.0; /* the Ensight data format */
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    for (int s = 0; s < enSet->scaN; ++s) {
        snprintf(enSet->fname, sizeof(EnStr), "%s.%s", enSet->bname, enSet->sca[s]);
        fp = Fopen(enSet->fname, "rb");
        Fread(enSet->str, sizeof(EnStr), 1, fp);
        for (int p = enSet->part[MIN], pnum = 1; p < enSet->part[MAX]; ++p, ++pnum) {
            Fread(enSet->str, sizeof(EnStr), 1, fp);
            Fread(&pnum, sizeof(int), 1, fp);
            Fread(enSet->str, sizeof(EnStr), 1, fp);
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
                        if (!InPartBox(k, j, i, part->ns[p])) {
                            continue;
                        }
                        /* data field initializer */
                        U = node[idx].U[TO];
                        Fread(&data, sizeof(EnReal), 1, fp);
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
        }
        fclose(fp);
    }
    return;
}
void ReadPolyDataEnsight(const Time *time, Geometry *const geo)
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
    EnSet enSet = { /* initialize environment */
        .rname = "geo_sph",
        .bname = {'\0'},
        .fname = {'\0'},
        .str = {'\0'},
        .fmt = "%s%05d",
        .gtag = "*****",
        .vtag = "*****",
        .dtype = "coordinates",
        .part = {0, 1},
        .scaN = 0,
        .sca = {{'\0'}},
        .vecN = 0,
        .vec = {{'\0'}},
    };
    snprintf(enSet.bname, sizeof(EnStr), enSet.fmt, enSet.rname, time->dataC);
    ReadPolyState(0, geo->sphN, geo, &enSet);
    return;
}
static void PolygonPolyDataReader(const Time *time, Geometry *const geo)
{
    EnSet enSet = { /* initialize environment */
        .rname = "geo_stl",
        .bname = {'\0'},
        .fname = {'\0'},
        .str = {'\0'},
        .fmt = "%s%05d",
        .gtag = "*****",
        .vtag = "*****",
        .dtype = "coordinates",
        .part = {geo->sphN, geo->totN},
        .scaN = 0,
        .sca = {{'\0'}},
        .vecN = 0,
        .vec = {{'\0'}},
    };
    snprintf(enSet.bname, sizeof(EnStr), enSet.fmt, enSet.rname, time->dataC);
    ReadPolygonPolyData(geo->sphN, geo->totN, geo, &enSet);
    return;
}
static void ReadPolygonPolyData(const int pm, const int pn, Geometry *const geo, EnSet *enSet)
{
    snprintf(enSet->fname, sizeof(EnStr), "%s.geo", enSet->bname);
    FILE *fp = Fopen(enSet->fname, "rb");
    EnReal data = 0.0; /* the Ensight data format */
    Polyhedron *poly = NULL;
    int ne = 0; /* total number of nodes in a part */
    Fread(enSet->str, sizeof(EnStr), 1, fp);
    Fread(enSet->str, sizeof(EnStr), 1, fp);
    Fread(enSet->str, sizeof(EnStr), 1, fp);
    Fread(enSet->str, sizeof(EnStr), 1, fp);
    Fread(enSet->str, sizeof(EnStr), 1, fp);
    for (int p = enSet->part[MIN], pnum = 1; p < enSet->part[MAX]; ++p, ++pnum) {
        poly = geo->poly + p;
        Fread(enSet->str, sizeof(EnStr), 1, fp);
        Fread(&pnum, sizeof(int), 1, fp);
        Fread(enSet->str, sizeof(EnStr), 1, fp);
        Sscanf(enSet->str, 3, "%d %d %d", &(poly->vertN), &(poly->edgeN), &(poly->faceN));
        Fread(enSet->str, sizeof(EnStr), 1, fp);
        Fread(&ne, sizeof(int), 1, fp);
        AllocatePolyhedronMemory(poly->vertN, poly->edgeN, poly->faceN, poly);
        poly->edgeN = 0; /* reset edge count before applying edge adding */
        for (int s = 0; s < DIMS; ++s) {
            for (int n = 0; n < poly->vertN; ++n) {
                Fread(&data, sizeof(EnReal), 1, fp);
                poly->v[n][s] = data;
            }
        }
        Fread(enSet->str, sizeof(EnStr), 1, fp);
        Fread(&ne, sizeof(int), 1, fp);
        for (int n = 0, m = 0; n < poly->faceN; ++n) {
            for (int s = 0; s < POLYN; ++s) {
                Fread(&m, sizeof(int), 1, fp);
                poly->f[n][s] = m - 1;
            }
            AddEdge(poly->f[n][0], poly->f[n][1], n, poly);
            AddEdge(poly->f[n][1], poly->f[n][2], n, poly);
            AddEdge(poly->f[n][2], poly->f[n][0], n, poly);
        }
        QuickSortEdge(poly->edgeN, poly->e);
    }
    ReadPolyState(pm, pn, geo, enSet);
    return;
}
static void ReadPolyState(const int pm, const int pn, Geometry *const geo, EnSet *enSet)
{
    snprintf(enSet->fname, sizeof(EnStr), "%s.state", enSet->bname);
    FILE *fp = Fopen(enSet->fname, "r");
    ReadPolyStateData(pm, pn, fp, geo);
    fclose(fp);
    return;
}
/* a good practice: end file with a newline */

