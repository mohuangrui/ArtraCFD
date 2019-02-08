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
#include "initialization.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "calculator.h"
#include "computational_geometry.h"
#include "immersed_boundary.h"
#include "boundary_treatment.h"
#include "data_stream.h"
#include "stl.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void InitializeSpaceData(Space *, const Model *);
static void InitializeFieldData(Space *, const Model *);
static void ApplyInitializer(const int, const Real [restrict],
        Real [restrict], const Partition *const, const Model *);
static void InitializeGeometryData(Geometry *const);
static void WritePolyMassProperty(const Geometry *const);
static void IdentifyGeometryState(Geometry *const);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void InitializeComputeDomain(Time *time, Space *space, const Model *model)
{
    if (0 == time->restart) { /* non restart */
        InitializeSpaceData(space, model);
    } else {
        ReadData(PROSD, time, space, model);
    }
    ComputeGeometryParameters(space->part.collapse, &(space->geo));
    WritePolyMassProperty(&(space->geo));
    ComputeGeometricField(space, model);
    TreatBoundary(TO, space, model);
    IdentifyGeometryState(&(space->geo));
    if (0 == time->restart) { /* non restart */
        WriteData(PROPT, time, space, model);
        WriteData(PROFC, time, space, model);
        WriteData(PROSD, time, space, model);
    }
    return;
}
static void InitializeSpaceData(Space *space, const Model *model)
{
    InitializeFieldData(space, model);
    InitializeGeometryData(&(space->geo));
    return;
}
/*
 * Initialize quantities for the entire domain
 * Exterior domains are initialized to unphysical values to avoid hiding
 * mistakes in boundary treatment and producing floating point exceptions.
 */
static void InitializeFieldData(Space *space, const Model *model)
{
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    RealVec pc = {0.0}; /* coordinates of current node */
    int idx = 0; /* linear array index math variable */
    for (int k = part->ns[PAL][Z][MIN]; k < part->ns[PAL][Z][MAX]; ++k) {
        for (int j = part->ns[PAL][Y][MIN]; j < part->ns[PAL][Y][MAX]; ++j) {
            for (int i = part->ns[PAL][X][MIN]; i < part->ns[PAL][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                node[idx].did = NONE;
                node[idx].fid = NONE;
                node[idx].lid = NONE;
                node[idx].gst = NONE;
                memset(node[idx].U, 1, DIMT * sizeof(*node[idx].U));
                if (!InPartBox(k, j, i, part->ns[PIN])) {
                    continue;
                }
                /* geometric field initializer */
                node[idx].did = 0;
                node[idx].fid = 0;
                node[idx].lid = 0;
                node[idx].gst = 0;
                /* data field initializer */
                pc[X] = MapPoint(i, part->domain[X][MIN], part->d[X], part->ng[X]);
                pc[Y] = MapPoint(j, part->domain[Y][MIN], part->d[Y], part->ng[Y]);
                pc[Z] = MapPoint(k, part->domain[Z][MIN], part->d[Z], part->ng[Z]);
                for (int n = 0; n < part->nIC; ++n) {
                    ApplyInitializer(n, pc, node[idx].U[TO], part, model);
                }
            }
        }
    }
    return;
}
static void ApplyInitializer(const int n, const Real pc[restrict], Real U[restrict],
        const Partition *const part, const Model *model)
{
    const Real zero = 0.0;
    const RealVec p1 = {part->posIC[n][0], part->posIC[n][1], part->posIC[n][2]};
    const RealVec p2 = {part->posIC[n][3], part->posIC[n][4], part->posIC[n][5]};
    const Real r = part->posIC[n][6];
    CalcVar var = {.t = zero, .x = pc[X], .y = pc[Y], .z = pc[Z], .ans = zero, .pi = PI};
    const Real Uo[DIMUo] = {
        ComputeExpression(&var, part->varIC[n][0]),
        ComputeExpression(&var, part->varIC[n][1]),
        ComputeExpression(&var, part->varIC[n][2]),
        ComputeExpression(&var, part->varIC[n][3]),
        ComputeExpression(&var, part->varIC[n][4])};
    const RealVec P1P2 = {p2[X] - p1[X], p2[Y] - p1[Y], p2[Z] - p1[Z]};
    const Real l2_P1P2 = Dot(P1P2, P1P2);
    RealVec P1Pc = {pc[X] - p1[X], pc[Y] - p1[Y], pc[Z] - p1[Z]};
    Real proj = zero; /* projection length */
    /* apply initial values for nodes that meets condition */
    int flag = 0; /* control flag for whether current node in the region */
    switch (part->typeIC[n]) {
        case ICGLOBAL:
            flag = 1;
            break;
        case ICPLANE:
            if (zero <= Dot(P1Pc, p2)) { /* on the normal direction or the plane */
                flag = 1;
            }
            break;
        case ICSPHERE:
            if (r * r >= Dot(P1Pc, P1Pc)) { /* in or on the sphere */
                flag = 1;
            }
            break;
        case ICBOX:
            P1Pc[X] = P1Pc[X] * (pc[X] - p2[X]);
            P1Pc[Y] = P1Pc[Y] * (pc[Y] - p2[Y]);
            P1Pc[Z] = P1Pc[Z] * (pc[Z] - p2[Z]);
            if ((zero >= P1Pc[X]) && (zero >= P1Pc[Y]) && (zero >= P1Pc[Z])) { /* in or on the box */
                flag = 1;
            }
            break;
        case ICCYLINDER:
            proj = Dot(P1Pc, P1P2);
            if ((zero > proj) || (l2_P1P2 < proj)) { /* outside the two ends */
                break;
            }
            proj = Dot(P1Pc, P1Pc) - proj * proj / l2_P1P2;
            if (r * r >= proj) { /* in or on the cylinder */
                flag = 1;
            }
            break;
        default:
            break;
    }
    if (1 == flag) { /* current node meets the condition */
        MapConservative(model->gamma, Uo, U);
    }
    return;
}
static void InitializeGeometryData(Geometry *const geo)
{
    FILE *fp = Fopen("artracfd.geo", "r");
    const char *fmtI = ParseFormat("%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg");
    /* read and process file line by line */
    String str = {'\0'}; /* store the current read line */
    String fname = {'\0'}; /* store the file name */
    while (NULL != fgets(str, sizeof str, fp)) {
        ParseCommand(str);
        if (0 == strncmp(str, "sphere state begin", sizeof str)) {
            ReadPolyStateData(0, geo->sphN, fp, geo);
            continue;
        }
        if (0 == strncmp(str, "polyhedron geometry begin", sizeof str)) {
            for (int n = geo->sphN; n < geo->totN; ++n) {
                Sread(fp, 1, "%s", fname);
                ReadStlFile(fname, geo->poly + n);
                ConvertPolyhedron(geo->poly + n);
            }
            continue;
        }
        if (0 == strncmp(str, "polyhedron state begin", sizeof str)) {
            ReadPolyStateData(geo->sphN, geo->totN, fp, geo);
            continue;
        }
        if (0 == strncmp(str, "polyhedron transform begin", sizeof str)) {
            const Real one = 1.0;
            const Real zero = 0.0;
            RealVec scale = {zero};
            RealVec angle = {zero};
            RealVec offset = {zero};
            for (int n = geo->sphN; n < geo->totN; ++n) {
                Sread(fp, 9, fmtI, scale + X, scale + Y, scale + Z,
                        angle + X, angle + Y, angle + Z, offset + X, offset + Y, offset + Z);
                if ((one == scale[X]) && (one == scale[Y]) && (one == scale[Z]) &&
                        (zero == angle[X]) && (zero == angle[Y]) && (zero == angle[Z]) &&
                        (zero == offset[X]) && (zero == offset[Y]) && (zero == offset[Z])) {
                    continue;
                }
                TransformPolyhedron(geo->poly[n].O, scale, angle, offset, geo->poly + n);
            }
            continue;
        }
    }
    fclose(fp);
    return;
}
static void WritePolyMassProperty(const Geometry *const geo)
{
    FILE *fp = Fopen("geo_mass_property.csv", "w");
    const char *fmtI = "  %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n";
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "#                          Polyhedron Mass Property                           -\n");
    fprintf(fp, "#                                                                             -\n");
    fprintf(fp, "# O[X,Y,Z], I[XX,YY,ZZ,XY,YZ,ZX], Box[X,Y,Z][MIN,MAX], rho, area, volume      -\n");
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    const Polyhedron *poly = NULL;
    for (int n = geo->sphN; n < geo->totN; ++n) {
        poly = geo->poly + n;
        fprintf(fp, fmtI,
                poly->O[X], poly->O[Y], poly->O[Z],
                poly->I[X][X], poly->I[Y][Y], poly->I[Z][Z],
                poly->I[X][Y], poly->I[Y][Z], poly->I[Z][X],
                poly->box[X][MIN], poly->box[Y][MIN], poly->box[Z][MIN],
                poly->box[X][MAX], poly->box[Y][MAX], poly->box[Z][MAX],
                poly->rho, poly->area, poly->volume);
    }
    fprintf(fp, "#------------------------------------------------------------------------------\n");
    fclose(fp);
    return;
}
static void IdentifyGeometryState(Geometry *const geo)
{
    Polyhedron *poly = NULL;
    const Real zero = 0.0;
    const Real rhoNoForce = 1.0e10;
    const Real rhoNoMove = 1.0e36;
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        if (rhoNoForce < poly->rho) { /* ignore surface force integration */
            poly->state = 2;
            if ((zero == poly->V[TO][X]) && (zero == poly->V[TO][Y]) && (zero == poly->V[TO][Z]) &&
                    (zero == poly->W[TO][X]) && (zero == poly->W[TO][Y]) && (zero == poly->W[TO][Z]) &&
                    (zero == poly->at[TN][X]) && (zero == poly->at[TN][Y]) && (zero == poly->at[TN][Z]) &&
                    (zero == poly->g[X]) && (zero == poly->g[Y]) && (zero == poly->g[Z]) &&
                    (zero == poly->ar[TN][X]) && (zero == poly->ar[TN][Y]) && (zero == poly->ar[TN][Z]) &&
                    (rhoNoMove < poly->rho)) { /* stationary geometry */
                poly->state = 1;
            }
        }
    }
    return;
}
/* a good practice: end file with a newline */

