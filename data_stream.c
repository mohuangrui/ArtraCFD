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
#include "data_stream.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include <float.h> /* size of floating point values */
#include "paraview.h"
#include "ensight.h"
#include "data_probe.h"
#include "commons.h"
/****************************************************************************
 * Function Pointers
 ****************************************************************************/
typedef void (*UnifiedDataWriter)(const Time *, const Space *, const Model *);
typedef void (*UnifiedDataReader)(Time *, Space *, const Model *);
typedef void (*StructuredDataWriter)(const Time *, const Space *, const Model *);
typedef void (*StructuredDataReader)(Time *, Space *, const Model *);
typedef void (*PolyDataWriter)(const Time *, const Geometry *const);
typedef void (*PolyDataReader)(const Time *, Geometry *const);
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void WriteSpaceData(const Time *, const Space *, const Model *);
static void ReadSpaceData(Time *, Space *, const Model *);
static void WriteFieldData(const Time *, const Space *, const Model *);
static void ReadFieldData(Time *, Space *, const Model *);
static void WriteGeometryData(const Time *, const Geometry *const);
static void ReadGeometryData(const Time *, Geometry *const);
static void WriteStateData(const Time *);
/****************************************************************************
 * Global Variables Definition with Private Scope
 ****************************************************************************/
static UnifiedDataWriter UnifiedWriteData[NPROBE] = {
    WritePointProbeData,
    WriteLineProbeData,
    WriteCurveProbeData,
    WriteSurfaceForceData,
    WriteSpaceData};
static UnifiedDataReader UnifiedReadData[NPROBE] = {
    ReadSpaceData,
    ReadSpaceData,
    ReadSpaceData,
    ReadSpaceData,
    ReadSpaceData};
static StructuredDataWriter WriteStructuredData[2] = {
    WriteStructuredDataParaview,
    WriteStructuredDataEnsight};
static StructuredDataReader ReadStructuredData[2] = {
    ReadStructuredDataParaview,
    ReadStructuredDataEnsight};
static PolyDataWriter WritePolyData[2] = {
    WritePolyDataParaview,
    WritePolyDataEnsight};
static PolyDataReader ReadPolyData[2] = {
    ReadPolyDataParaview,
    ReadPolyDataEnsight};
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void WriteData(const int n, const Time *time, const Space *space, const Model *model)
{
    UnifiedWriteData[n](time, space, model);
    return;
}
void ReadData(const int n, Time *time, Space *space, const Model *model)
{
    UnifiedReadData[n](time, space, model);
    return;
}
static void WriteSpaceData(const Time *time, const Space *space, const Model *model)
{
    WriteFieldData(time, space, model);
    WriteGeometryData(time, &(space->geo));
    WriteStateData(time);
    return;
}
static void ReadSpaceData(Time *time, Space *space, const Model *model)
{
    ReadFieldData(time, space, model);
    ReadGeometryData(time, &(space->geo));
    return;
}
static void WriteFieldData(const Time *time, const Space *space, const Model *model)
{
    WriteStructuredData[time->dataStreamer](time, space, model);
    return;
}
static void ReadFieldData(Time *time, Space *space, const Model *model)
{
    ReadStructuredData[time->dataStreamer](time, space, model);
    return;
}
static void WriteGeometryData(const Time *time, const Geometry *const geo)
{
    if (0 == geo->totN) {
        return;
    }
    WritePolyData[time->dataStreamer](time, geo);
    return;
}
static void ReadGeometryData(const Time *time, Geometry *const geo)
{
    if (0 == geo->totN) {
        return;
    }
    ReadPolyData[time->dataStreamer](time, geo);
    return;
}
static void WriteStateData(const Time *time)
{
    const char *fname = "artracfd.log";
    FILE *fp = Fopen(fname, "w");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "%d  # completion checkpoint\n", (time->now == time->end) || (time->stepC == time->stepN));
    fprintf(fp, "%d  # data checkpoint\n", time->dataC);
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fclose(fp);
    return;
}
void WritePolyStateData(const int pm, const int pn, FILE *fp, const Geometry *const geo)
{
    const char *fmtI = "  %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %d\n";
    const char *fmtJ = "  %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n";
    const Polyhedron *poly = NULL;
    for (int n = pm; n < pn; ++n) {
        poly = geo->poly + n;
        fprintf(fp, fmtI,
                poly->O[X], poly->O[Y], poly->O[Z], poly->r,
                poly->V[TO][X], poly->V[TO][Y], poly->V[TO][Z],
                poly->W[TO][X], poly->W[TO][Y], poly->W[TO][Z],
                poly->rho, poly->T, poly->cf,
                poly->area, poly->volume, poly->mid);
        fprintf(fp, fmtJ,
                poly->at[TO][X], poly->at[TO][Y], poly->at[TO][Z],
                poly->ar[TO][X], poly->ar[TO][Y], poly->ar[TO][Z],
                poly->at[TN][X], poly->at[TN][Y], poly->at[TN][Z],
                poly->g[X], poly->g[Y], poly->g[Z],
                poly->ar[TN][X], poly->ar[TN][Y], poly->ar[TN][Z],
                poly->to);
    }
    return;
}
void ReadPolyStateData(const int pm, const int pn, FILE *fp, Geometry *const geo)
{
    const char *fmtI = ParseFormat("%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %d");
    const char *fmtJ = ParseFormat("%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg");
    Polyhedron *poly  = NULL;
    const Real zero = 0.0;
    for (int n = pm; n < pn; ++n) {
        poly = geo->poly + n;
        Sread(fp, 16, fmtI,
                &(poly->O[X]), &(poly->O[Y]), &(poly->O[Z]), &(poly->r),
                &(poly->V[TO][X]), &(poly->V[TO][Y]), &(poly->V[TO][Z]),
                &(poly->W[TO][X]), &(poly->W[TO][Y]), &(poly->W[TO][Z]),
                &(poly->rho), &(poly->T), &(poly->cf),
                &(poly->area), &(poly->volume), &(poly->mid));
        Sread(fp, 16, fmtJ,
                &(poly->at[TO][X]), &(poly->at[TO][Y]), &(poly->at[TO][Z]),
                &(poly->ar[TO][X]), &(poly->ar[TO][Y]), &(poly->ar[TO][Z]),
                &(poly->at[TN][X]), &(poly->at[TN][Y]), &(poly->at[TN][Z]),
                &(poly->g[X]), &(poly->g[Y]), &(poly->g[Z]),
                &(poly->ar[TN][X]), &(poly->ar[TN][Y]), &(poly->ar[TN][Z]),
                &(poly->to));
        if (zero >= poly->to) {
            poly->to = FLT_MAX;
        }
        if (geo->sphN > n) {
            poly->faceN = 0; /* analytical polyhedron tag */
            poly->facet = NULL;
        }
        for (int s = 0; s < DIMS; ++s) {
            poly->V[TN][s] = poly->V[TO][s];
            poly->W[TN][s] = poly->W[TO][s];
        }
    }
    return;
}
/* a good practice: end file with a newline */

