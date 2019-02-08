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
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void InitializeTransientCaseFile(EnSet *);
static void WriteCaseFile(const Time *, EnSet *);
static void WriteGeometryFile(const Space *, EnSet *);
static void WriteStructuredData(const Space *, const Model *, EnSet *);
static void PointPolyDataWriter(const Time *, const Geometry *const);
static void WritePointPolyData(const int, const int, const Geometry *const, EnSet *);
static void PolygonPolyDataWriter(const Time *, const Geometry *const);
static void WritePolygonPolyData(const int, const int, const Geometry *const, EnSet *);
static void WritePolyVariable(const int, const int, const Geometry *const, EnSet *);
static void WritePolyState(const int, const int, const Geometry *const, EnSet *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void WriteStructuredDataEnsight(const Time *time, const Space *space, const Model *model)
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
        .scaN = 7,
        .sca = {"rho", "u", "v", "w", "p", "T", "did"},
        .vecN = 1,
        .vec = {"Vel"},
    };
    snprintf(enSet.bname, sizeof(EnStr), enSet.fmt, enSet.rname, time->dataC);
    if (0 == time->stepC) { /* initialization step */
        InitializeTransientCaseFile(&enSet);
        WriteGeometryFile(space, &enSet);
    }
    WriteCaseFile(time, &enSet);
    WriteStructuredData(space, model, &enSet);
    return;
}
static void InitializeTransientCaseFile(EnSet *enSet)
{
    snprintf(enSet->fname, sizeof(EnStr), "%s.case", enSet->rname);
    FILE *fp = Fopen(enSet->fname, "w");
    fprintf(fp, "FORMAT\n");
    fprintf(fp, "type: ensight gold\n");
    fprintf(fp, "\n");
    fprintf(fp, "GEOMETRY\n");
    fprintf(fp, "model: %s%s.geo\n", enSet->rname, enSet->gtag);
    fprintf(fp, "\n");
    fprintf(fp, "VARIABLE\n");
    for (int n = 0; n < enSet->scaN; ++n) {
        fprintf(fp, "scalar per node:  1  %3s  %s%s.%s\n",
                enSet->sca[n], enSet->rname, enSet->vtag, enSet->sca[n]);
    }
    for (int n = 0; n < enSet->vecN; ++n) {
        fprintf(fp, "vector per node:  1  %3s  %s%s.%s\n",
                enSet->vec[n], enSet->rname, enSet->vtag, enSet->vec[n]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "TIME\n");
    fprintf(fp, "time set: 1\n");
    fprintf(fp, "number of steps:          0          \n");
    fprintf(fp, "filename start number:    0\n");
    fprintf(fp, "filename increment:       1\n");
    fprintf(fp, "time values:  ");
    fclose(fp);
    return;
}
static void WriteCaseFile(const Time *time, EnSet *enSet)
{
    snprintf(enSet->fname, sizeof(EnStr), "%s.case", enSet->bname);
    FILE *fp = Fopen(enSet->fname, "w");
    fprintf(fp, "FORMAT\n");
    fprintf(fp, "type: ensight gold\n");
    fprintf(fp, "\n");
    fprintf(fp, "GEOMETRY\n");
    if ('\0' == *enSet->gtag) {
        fprintf(fp, "model: %s.geo\n", enSet->rname);
    } else {
        fprintf(fp, "model: %s.geo\n", enSet->bname);
    }
    fprintf(fp, "\n");
    fprintf(fp, "VARIABLE\n");
    fprintf(fp, "constant per case:  Time  %.6g\n", time->now);
    fprintf(fp, "constant per case:  Step  %d\n", time->stepC);
    for (int n = 0; n < enSet->scaN; ++n) {
        fprintf(fp, "scalar per node:     %3s  %s.%s\n",
                enSet->sca[n], enSet->bname, enSet->sca[n]);
    }
    for (int n = 0; n < enSet->vecN; ++n) {
        fprintf(fp, "vector per node:     %3s  %s.%s\n",
                enSet->vec[n], enSet->bname, enSet->vec[n]);
    }
    fprintf(fp, "\n");
    fclose(fp);
    /* add case to the transient case file */
    snprintf(enSet->fname, sizeof(EnStr), "%s.case", enSet->rname);
    fp = Fopen(enSet->fname, "r+");
    /* seek the target line for adding information */
    ReadInLine(fp, "time set: 1");
    fprintf(fp, "number of steps:          %d", (time->dataC + 1));
    /* add the time flag of current export to the transient case */
    fseek(fp, 0, SEEK_END); /* seek to the end of file */
    if ((time->dataC % 5) == 0) { /* print to a new line every x outputs */
        fprintf(fp, "\n");
    }
    fprintf(fp, "%.6g ", time->now);
    fclose(fp);
    return;
}
static void WriteGeometryFile(const Space *space, EnSet *enSet)
{
    /*
     * Write the geometry file in Binary Form.
     * Maximums: maximum number of nodes in a part is 2GB.
     */
    snprintf(enSet->fname, sizeof(EnStr), "%s.geo", enSet->rname);
    FILE *fp = Fopen(enSet->fname, "wb");
    EnReal data = 0.0; /* the Ensight data format */
    const Partition *const part = &(space->part);
    IntVec ne = {0}; /* i, j, k node number in each part */
    /* description at the beginning */
    strncpy(enSet->str, "C Binary", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    strncpy(enSet->str, "Ensight Geometry File", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    strncpy(enSet->str, "Written by ArtraCFD", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    /* node id and extents settings */
    strncpy(enSet->str, "node id off", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    strncpy(enSet->str, "element id off", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    for (int p = enSet->part[MIN], pnum = 1; p < enSet->part[MAX]; ++p, ++pnum) {
        strncpy(enSet->str, "part", sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        fwrite(&pnum, sizeof(int), 1, fp);
        snprintf(enSet->str, sizeof(EnStr), "part %d", p);
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        strncpy(enSet->str, enSet->dtype, sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        ne[X] = part->ns[p][X][MAX] - part->ns[p][X][MIN];
        ne[Y] = part->ns[p][Y][MAX] - part->ns[p][Y][MIN];
        ne[Z] = part->ns[p][Z][MAX] - part->ns[p][Z][MIN];
        fwrite(ne, sizeof(int), 3, fp);
        for (int s = 0; s < DIMS; ++s) {
            for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
                for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
                    for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                        ne[X] = i; ne[Y] = j; ne[Z] = k;
                        data = MapPoint(ne[s], part->domain[s][MIN], part->d[s], part->ng[s]);
                        fwrite(&data, sizeof(EnReal), 1, fp);
                    }
                }
            }
        }
    }
    fclose(fp);
    return;
}
/*
 * The values for each node of the structured block are output in
 * the same IJK order as the coordinates. (The number of nodes in the
 * part are obtained from the corresponding geometry file.)
 */
static void WriteStructuredData(const Space *space, const Model *model, EnSet *enSet)
{
    FILE *fp = NULL;
    EnReal data = 0.0; /* the Ensight data format */
    const Partition *const part = &(space->part);
    const Node *const node = space->node;
    const Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    for (int s = 0; s < enSet->scaN; ++s) {
        snprintf(enSet->fname, sizeof(EnStr), "%s.%s", enSet->bname, enSet->sca[s]);
        fp = Fopen(enSet->fname, "wb");
        /* first line description per file */
        strncpy(enSet->str, "scalar variable", sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        for (int p = enSet->part[MIN], pnum = 1; p < enSet->part[MAX]; ++p, ++pnum) {
            /* binary file format */
            strncpy(enSet->str, "part", sizeof(EnStr));
            fwrite(enSet->str, sizeof(EnStr), 1, fp);
            fwrite(&pnum, sizeof(int), 1, fp);
            strncpy(enSet->str, enSet->dtype, sizeof(EnStr));
            fwrite(enSet->str, sizeof(EnStr), 1, fp);
            /* now output the scalar value at each node in current part */
            for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
                for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
                    for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
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
                            default:
                                break;
                        }
                        fwrite(&data, sizeof(EnReal), 1, fp);
                    }
                }
            }
        }
        fclose(fp);
    }
    for (int s = 0; s < enSet->vecN; ++s) {
        snprintf(enSet->fname, sizeof(EnStr), "%s.%s", enSet->bname, enSet->vec[s]);
        fp = Fopen(enSet->fname, "wb");
        /* binary file format */
        strncpy(enSet->str, "vector variable", sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        for (int p = enSet->part[MIN], pnum = 1; p < enSet->part[MAX]; ++p, ++pnum) {
            strncpy(enSet->str, "part", sizeof(EnStr));
            fwrite(enSet->str, sizeof(EnStr), 1, fp);
            fwrite(&pnum, sizeof(int), 1, fp);
            strncpy(enSet->str, enSet->dtype, sizeof(EnStr));
            fwrite(enSet->str, sizeof(EnStr), 1, fp);
            for (int n = 1; n < 4; ++n) {
                for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
                    for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
                        for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                            idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                            U = node[idx].U[TO];
                            data = U[n] / U[0];
                            fwrite(&data, sizeof(EnReal), 1, fp);
                        }
                    }
                }
            }
        }
        fclose(fp);
    }
    return;
}
void WritePolyDataEnsight(const Time *time, const Geometry *const geo)
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
        .scaN = 2,
        .sca = {"r", "did"},
        .vecN = 1,
        .vec = {"Vel"},
    };
    snprintf(enSet.bname, sizeof(EnStr), enSet.fmt, enSet.rname, time->dataC);
    if (0 == time->stepC) { /* initialization step */
        InitializeTransientCaseFile(&enSet);
    }
    WriteCaseFile(time, &enSet);
    WritePointPolyData(0, geo->sphN, geo, &enSet);
    return;
}
static void WritePointPolyData(const int pm, const int pn, const Geometry *const geo, EnSet *enSet)
{
    snprintf(enSet->fname, sizeof(EnStr), "%s.geo", enSet->bname);
    FILE *fp = Fopen(enSet->fname, "wb");
    EnReal data = 0.0; /* the Ensight data format */
    int ne = 0; /* total number of nodes in a part */
    /* description at the beginning */
    strncpy(enSet->str, "C Binary", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    strncpy(enSet->str, "Ensight Geometry File", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    strncpy(enSet->str, "Written by ArtraCFD", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    /* node id and extents settings */
    strncpy(enSet->str, "node id off", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    strncpy(enSet->str, "element id off", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    for (int p = enSet->part[MIN], pnum = 1; p < enSet->part[MAX]; ++p, ++pnum) {
        strncpy(enSet->str, "part", sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        fwrite(&pnum, sizeof(int), 1, fp);
        snprintf(enSet->str, sizeof(EnStr), "%d", pn - pm);
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        strncpy(enSet->str, enSet->dtype, sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        ne = pn - pm;
        fwrite(&ne, sizeof(int), 1, fp);
        for (int s = 0; s < DIMS; ++s) {
            for (int n = pm; n < pn; ++n) {
                data = geo->poly[n].O[s];
                fwrite(&data, sizeof(EnReal), 1, fp);
            }
        }
        strncpy(enSet->str, "point", sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        ne = pn - pm;
        fwrite(&ne, sizeof(int), 1, fp);
        for (int n = pm, m = 1; n < pn; ++n, ++m) {
            fwrite(&m, sizeof(int), 1, fp);
        }
    }
    fclose(fp);
    WritePolyVariable(pm, pn, geo, enSet);
    WritePolyState(pm, pn, geo, enSet);
    return;
}
static void PolygonPolyDataWriter(const Time *time, const Geometry *const geo)
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
    if (0 == time->stepC) { /* initialization step */
        InitializeTransientCaseFile(&enSet);
    }
    WriteCaseFile(time, &enSet);
    WritePolygonPolyData(geo->sphN, geo->totN, geo, &enSet);
    return;
}
static void WritePolygonPolyData(const int pm, const int pn, const Geometry *const geo, EnSet *enSet)
{
    snprintf(enSet->fname, sizeof(EnStr), "%s.geo", enSet->bname);
    FILE *fp = Fopen(enSet->fname, "wb");
    EnReal data = 0.0; /* the Ensight data format */
    const Polyhedron *poly = NULL;
    int ne = 0; /* total number of nodes in a part */
    /* description at the beginning */
    strncpy(enSet->str, "C Binary", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    strncpy(enSet->str, "Ensight Geometry File", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    strncpy(enSet->str, "Written by ArtraCFD", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    /* node id and extents settings */
    strncpy(enSet->str, "node id off", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    strncpy(enSet->str, "element id off", sizeof(EnStr));
    fwrite(enSet->str, sizeof(EnStr), 1, fp);
    for (int p = enSet->part[MIN], pnum = 1; p < enSet->part[MAX]; ++p, ++pnum) {
        poly = geo->poly + p;
        strncpy(enSet->str, "part", sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        fwrite(&pnum, sizeof(int), 1, fp);
        snprintf(enSet->str, sizeof(EnStr), "%d %d %d",
                poly->vertN, poly->edgeN, poly->faceN);
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        strncpy(enSet->str, enSet->dtype, sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        ne = poly->vertN;
        fwrite(&ne, sizeof(int), 1, fp);
        for (int s = 0; s < DIMS; ++s) {
            for (int n = 0; n < poly->vertN; ++n) {
                data = poly->v[n][s];
                fwrite(&data, sizeof(EnReal), 1, fp);
            }
        }
        strncpy(enSet->str, "tria3", sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        ne = poly->faceN;
        fwrite(&ne, sizeof(int), 1, fp);
        for (int n = 0, m = 0; n < poly->faceN; ++n) {
            for (int s = 0; s < POLYN; ++s) {
                m = poly->f[n][s] + 1;
                fwrite(&m, sizeof(int), 1, fp);
            }
        }
    }
    fclose(fp);
    WritePolyState(pm, pn, geo, enSet);
    return;
}
static void WritePolyVariable(const int pm, const int pn, const Geometry *const geo, EnSet *enSet)
{
    FILE *fp = NULL;
    EnReal data = 0.0; /* the Ensight data format */
    for (int s = 0; s < enSet->scaN; ++s) {
        snprintf(enSet->fname, sizeof(EnStr), "%s.%s", enSet->bname, enSet->sca[s]);
        fp = Fopen(enSet->fname, "wb");
        /* first line description per file */
        strncpy(enSet->str, "scalar variable", sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        for (int p = enSet->part[MIN], pnum = 1; p < enSet->part[MAX]; ++p, ++pnum) {
            /* binary file format */
            strncpy(enSet->str, "part", sizeof(EnStr));
            fwrite(enSet->str, sizeof(EnStr), 1, fp);
            fwrite(&pnum, sizeof(int), 1, fp);
            strncpy(enSet->str, enSet->dtype, sizeof(EnStr));
            fwrite(enSet->str, sizeof(EnStr), 1, fp);
            /* now output the scalar value at each node in current part */
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
                fwrite(&data, sizeof(EnReal), 1, fp);
            }
        }
        fclose(fp);
    }
    for (int s = 0; s < enSet->vecN; ++s) {
        snprintf(enSet->fname, sizeof(EnStr), "%s.%s", enSet->bname, enSet->vec[s]);
        fp = Fopen(enSet->fname, "wb");
        /* binary file format */
        strncpy(enSet->str, "vector variable", sizeof(EnStr));
        fwrite(enSet->str, sizeof(EnStr), 1, fp);
        for (int p = enSet->part[MIN], pnum = 1; p < enSet->part[MAX]; ++p, ++pnum) {
            strncpy(enSet->str, "part", sizeof(EnStr));
            fwrite(enSet->str, sizeof(EnStr), 1, fp);
            fwrite(&pnum, sizeof(int), 1, fp);
            strncpy(enSet->str, enSet->dtype, sizeof(EnStr));
            fwrite(enSet->str, sizeof(EnStr), 1, fp);
            for (int n = 0; n < DIMS; ++n) {
                for (int m = pm; m < pn; ++m) {
                    data = geo->poly[m].V[TO][n];
                    fwrite(&data, sizeof(EnReal), 1, fp);
                }
            }
        }
        fclose(fp);
    }
    return;
}
static void WritePolyState(const int pm, const int pn, const Geometry *const geo, EnSet *enSet)
{
    snprintf(enSet->fname, sizeof(EnStr), "%s.state", enSet->bname);
    FILE *fp = Fopen(enSet->fname, "w");
    WritePolyStateData(pm, pn, fp, geo);
    fclose(fp);
    return;
}
/* a good practice: end file with a newline */

