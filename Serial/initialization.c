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
#include "computational_geometry.h"
#include "immersed_boundary.h"
#include "boundary_treatment.h"
#include "data_stream.h"
#include "data_probe.h"
#include "paraview.h"
#include "stl.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int GlobalInitialization(Space *);
static int InitializeFieldData(Space *, const Model *);
static int ApplyRegionalInitializer(const int, Space *, const Model *);
static int InitializeGeometryData(Geometry *);
static int WritePolyhedronMassProperty(const Geometry *);
static int IdentifyGeometryState(Geometry *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int InitializeComputationalDomain(Time *time, Space *space, const Model *model)
{
    GlobalInitialization(space);
    if (0 == time->restart) { /* non restart */
        InitializeFieldData(space, model);
        InitializeGeometryData(&(space->geo));
    } else {
        ReadFieldData(time, space, model);
        ReadGeometryData(time, &(space->geo));
    }
    ComputeGeometryParameters(space->part.collapse, &(space->geo));
    WritePolyhedronMassProperty(&(space->geo));
    ComputeGeometryDomain(space, model);
    BoundaryConditionsAndTreatments(TO, space, model);
    IdentifyGeometryState(&(space->geo));
    if (0 == time->restart) { /* non restart */
        WriteSurfaceForceData(time, space);
        WriteFieldData(time, space, model);
        WriteGeometryData(time, &(space->geo));
        WriteFieldDataAtPointProbes(time, space, model);
    }
    return 0;
}
/*
 * Initialize non zero quantities for the global domain.
 * After this global initialization, future computation only need
 * to focus on the interior computational domain. Physical quantities
 * are initialized into unphysical values to avoid hiding potential
 * mistakes in mesh generation and boundary treatment.
 */
static int GlobalInitialization(Space *space)
{
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    const int idxMax = part->n[X] * part->n[Y] * part->n[Z];
    for (int idx = 0; idx < idxMax; ++idx) {
        node[idx].gid = NONE;
        node[idx].fid = NONE;
        node[idx].lid = NONE;
        node[idx].gst = NONE;
    }
    return 0;
}
static int InitializeFieldData(Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    /* extract global initial values */
    const Real Uo[DIMUo] = {
        part->valueIC[0][ENTRYIC-5],
        part->valueIC[0][ENTRYIC-4],
        part->valueIC[0][ENTRYIC-3],
        part->valueIC[0][ENTRYIC-2],
        part->valueIC[0][ENTRYIC-1]};
    /*
     * Initialize the interior field
     */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                U = node[idx].U[TO];
                ConservativeByPrimitive(model->gamma, Uo, U);
            }
        }
    }
    /*
     * Regional initializer for specific regions
     */
    for (int n = 1; n < part->countIC; ++n) {
        ApplyRegionalInitializer(n, space, model);
    }
    return 0;
}
/*
 * The handling of regional initialization for specific region is achieved
 * through the cooperation of three data structures:
 * The countIC counts the number of initializers.
 * The typeIC array keeps a list of the types of regional initialization.
 * The valueIC array stored the information of the corresponding IC type.
 */
static int ApplyRegionalInitializer(const int n, Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    Real *restrict U = NULL;
    const Real zero = 0.0;
    int idx = 0; /* linear array index math variable */
    /*
     * Acquire the specialized information data entries
     */
    const RealVec p1 = {
        part->valueIC[n][0],
        part->valueIC[n][1],
        part->valueIC[n][2]};
    const RealVec p2 = {
        part->valueIC[n][3],
        part->valueIC[n][4],
        part->valueIC[n][5]};
    const Real r = part->valueIC[n][6];
    const Real Uo[DIMUo] = {
        part->valueIC[n][ENTRYIC-5],
        part->valueIC[n][ENTRYIC-4],
        part->valueIC[n][ENTRYIC-3],
        part->valueIC[n][ENTRYIC-2],
        part->valueIC[n][ENTRYIC-1]};
    const RealVec P1P2 = {
        p2[X] - p1[X],
        p2[Y] - p1[Y],
        p2[Z] - p1[Z]};
    const Real l2_P1P2 = Dot(P1P2, P1P2);
    /*
     * Apply initial values for nodes that meets condition
     */
    int flag = 0; /* control flag for whether current node in the region */
    RealVec pc = {zero}; /* coordinates of current node */
    RealVec P1Pc = {zero}; /* position vector */
    Real proj = zero; /* projection length */
    for (int k = part->ns[PIN][Z][MIN]; k < part->ns[PIN][Z][MAX]; ++k) {
        for (int j = part->ns[PIN][Y][MIN]; j < part->ns[PIN][Y][MAX]; ++j) {
            for (int i = part->ns[PIN][X][MIN]; i < part->ns[PIN][X][MAX]; ++i) {
                pc[X] = PointSpace(i, part->domain[X][MIN], part->d[X], part->ng);
                pc[Y] = PointSpace(j, part->domain[Y][MIN], part->d[Y], part->ng);
                pc[Z] = PointSpace(k, part->domain[Z][MIN], part->d[Z], part->ng);
                P1Pc[X] = pc[X] - p1[X];
                P1Pc[Y] = pc[Y] - p1[Y];
                P1Pc[Z] = pc[Z] - p1[Z];
                flag = 0; /* always initialize flag to zero */
                switch (part->typeIC[n]) {
                    case ICPLANE:
                        if (zero <= Dot(P1Pc, p2)) { /* on the normal direction or the plane */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case ICSPHERE:
                        if (r * r >= Dot(P1Pc, P1Pc)) { /* in or on the sphere */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case ICBOX:
                        P1Pc[X] = P1Pc[X] * (pc[X] - p2[X]);
                        P1Pc[Y] = P1Pc[Y] * (pc[Y] - p2[Y]);
                        P1Pc[Z] = P1Pc[Z] * (pc[Z] - p2[Z]);
                        if ((zero >= P1Pc[X]) && (zero >= P1Pc[Y]) && (zero >= P1Pc[Z])) { /* in or on the box */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case ICCYLINDER:
                        proj = Dot(P1Pc, P1P2);
                        if ((zero > proj) || (l2_P1P2 < proj)) { /* outside the two ends */
                            break;
                        }
                        proj = Dot(P1Pc, P1Pc) - proj * proj / l2_P1P2;
                        if (r * r >= proj) { /* in or on the cylinder */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    default:
                        break;
                }
                if (1 == flag) { /* current node meets the condition */
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    U = node[idx].U[TO];
                    ConservativeByPrimitive(model->gamma, Uo, U);
                }
            }
        }
    }
    return 0;
}
static int InitializeGeometryData(Geometry *geo)
{
    FILE *filePointer = fopen("artracfd.geo", "r");
    if (NULL == filePointer) {
        FatalError("failed to open file: artracfd.geo...");
    }
    char formatIX[45] = "%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatIX, "%g, %g, %g, %g, %g, %g, %g, %g, %g", sizeof formatIX); /* float type */
    }
    /* read and process file line by line */
    String currentLine = {'\0'}; /* store the current read line */
    int nscan = 0; /* read conversion count */
    String fileName = {'\0'}; /* store the file name */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "sphere state begin", sizeof currentLine)) {
            ReadPolyhedronStateData(0, geo->sphN, filePointer, geo);
            continue;
        }
        if (0 == strncmp(currentLine, "polyhedron geometry begin", sizeof currentLine)) {
            for (int n = geo->sphN; n < geo->totN; ++n) {
                Fgets(currentLine, sizeof currentLine, filePointer);
                nscan = sscanf(currentLine, "%s", fileName);
                VerifyReadConversion(nscan, 1);
                ReadStlFile(fileName, geo->poly + n);
                ConvertPolyhedron(geo->poly + n);
            }
            continue;
        }
        if (0 == strncmp(currentLine, "polyhedron state begin", sizeof currentLine)) {
            ReadPolyhedronStateData(geo->sphN, geo->totN, filePointer, geo);
            continue;
        }
        if (0 == strncmp(currentLine, "polyhedron transform begin", sizeof currentLine)) {
            const Real one = 1.0;
            const Real zero = 0.0;
            RealVec scale = {zero};
            RealVec angle = {zero};
            RealVec offset = {zero};
            for (int n = geo->sphN; n < geo->totN; ++n) {
                Fgets(currentLine, sizeof currentLine, filePointer);
                nscan = sscanf(currentLine, formatIX, scale + X, scale + Y, scale + Z, 
                        angle + X, angle + Y, angle + Z, offset + X, offset + Y, offset + Z);
                VerifyReadConversion(nscan, 9);
                if ((one == scale[X]) && (one == scale[Y]) && (one == scale[Z]) &&
                        (zero == angle[X]) && (zero == angle[Y]) && (zero == angle[Z]) &&
                        (zero == offset[X]) && (zero == offset[Y]) && (zero == offset[Z])) {
                    continue;
                }
                Transformation(geo->poly[n].O, scale, angle, offset, geo->poly + n);
            }
            continue;
        }
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int WritePolyhedronMassProperty(const Geometry *geo)
{
    FILE *filePointer = fopen("geo_mass_property.csv", "w");
    if (NULL == filePointer) {
        FatalError("failed to open data file...");
    }
    const char format[110] = "  %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n";
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#                          Polyhedron Mass Property                           -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "# O[X,Y,Z], I[XX,YY,ZZ,XY,YZ,ZX], Box[X,Y,Z][MIN,MAX], rho, area, volume      -\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    const Polyhedron *poly = NULL;
    for (int n = geo->sphN; n < geo->totN; ++n) {
        poly = geo->poly + n;
        fprintf(filePointer, format,
                poly->O[X], poly->O[Y], poly->O[Z], 
                poly->I[X][X], poly->I[Y][Y], poly->I[Z][Z],
                poly->I[X][Y], poly->I[Y][Z], poly->I[Z][X],
                poly->box[X][MIN], poly->box[Y][MIN], poly->box[Z][MIN],
                poly->box[X][MAX], poly->box[Y][MAX], poly->box[Z][MAX],
                poly->rho, poly->area, poly->volume);
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * Identify geometry state
 *
 * Identify some special geometry state to simplfy computation.
 * Should be done after the initial computation of geometry domain.
 */
static int IdentifyGeometryState(Geometry *geo)
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
    return 0;
}
/* a good practice: end file with a newline */

