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
#include "data_probe.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* support for abs operation */
#include "cfd_commons.h"
#include "immersed_boundary.h"
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void WriteFieldDataAtPointProbes(const Time *time, const Space *space, const Model *model)
{
    if (0 == time->pointProbeN) {
        return;
    }
    FILE *filePointer = NULL;
    String fileName = {'\0'};
    const Partition *restrict part = &(space->part);
    const Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0};
    const IntVec nMin = {part->ng, part->ng, part->ng};
    const IntVec nMax = {part->n[X] - part->ng, part->n[Y] - part->ng, part->n[Z] - part->ng};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const int ng = part->ng;
    RealVec p1 = {0.0};
    int i = 0, j = 0, k = 0;
    for (int n = 0; n < time->pointProbeN; ++n) {
        snprintf(fileName, sizeof(fileName), "%s%03d.csv", "point_probe_", n + 1);
        filePointer = fopen(fileName, "a");
        if (NULL == filePointer) {
            FatalError("failed to write data at probes...");
        }
        if (0 == time->stepC) { /* this is the initialization step */
            fprintf(filePointer, "# time, rho, u, v, w, p, T\n"); 
        }
        p1[X] = time->pp[n][0];
        p1[Y] = time->pp[n][1];
        p1[Z] = time->pp[n][2];
        i = ValidNodeSpace(NodeSpace(p1[X], sMin[X], dd[X], ng), nMin[X], nMax[X]);
        j = ValidNodeSpace(NodeSpace(p1[Y], sMin[Y], dd[Y], ng), nMin[Y], nMax[Y]);
        k = ValidNodeSpace(NodeSpace(p1[Z], sMin[Z], dd[Z], ng), nMin[Z], nMax[Z]);
        idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
        PrimitiveByConservative(model->gamma, model->gasR, node[idx].U[TO], Uo);
        fprintf(filePointer, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                time->now, Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]); 
        fclose(filePointer); /* close current opened file */
    }
    return;
}
void WriteFieldDataAtLineProbes(const Time *time, const Space *space, const Model *model)
{
    if (0 == time->lineProbeN) {
        return;
    }
    FILE *filePointer = NULL;
    String fileName = {'\0'};
    const Partition *restrict part = &(space->part);
    const Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    int idxOld = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0};
    const IntVec nMin = {part->ng, part->ng, part->ng};
    const IntVec nMax = {part->n[X] - part->ng, part->n[Y] - part->ng, part->n[Z] - part->ng};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const int ng = part->ng;
    RealVec p1 = {0.0};
    RealVec p2 = {0.0};
    RealVec dl = {0.0};
    int stepN = 0;
    int i = 0, j = 0, k = 0;
    for (int n = 0; n < time->lineProbeN; ++n) {
        snprintf(fileName, sizeof(fileName), "%s%03d_%05d.csv", "line_probe_", n + 1, time->stepC);
        filePointer = fopen(fileName, "w");
        if (NULL == filePointer) {
            FatalError("failed to write data at probes...");
        }
        fprintf(filePointer, "# x, y, z, rho, u, v, w, p, T <time=%.6g>\n", time->now); 
        p1[X] = time->lp[n][0];
        p1[Y] = time->lp[n][1];
        p1[Z] = time->lp[n][2];
        p2[X] = time->lp[n][3];
        p2[Y] = time->lp[n][4];
        p2[Z] = time->lp[n][5];
        stepN = MaxInt(time->lp[n][6] - 1, 1);
        dl[X] = (p2[X] - p1[X]) / (Real)(stepN);
        dl[Y] = (p2[Y] - p1[Y]) / (Real)(stepN);
        dl[Z] = (p2[Z] - p1[Z]) / (Real)(stepN);
        idxOld = -1; /* used to avoid repeating node for tiny step sizes */
        for (int m = 0; m <= stepN; ++m) {
            i = ValidNodeSpace(NodeSpace(p1[X] + m * dl[X], sMin[X], dd[X], ng), nMin[X], nMax[X]);
            j = ValidNodeSpace(NodeSpace(p1[Y] + m * dl[Y], sMin[Y], dd[Y], ng), nMin[Y], nMax[Y]);
            k = ValidNodeSpace(NodeSpace(p1[Z] + m * dl[Z], sMin[Z], dd[Z], ng), nMin[Z], nMax[Z]);
            idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
            if (idxOld == idx) {
                continue;
            }
            idxOld = idx; /* record */
            p2[X] = PointSpace(i, sMin[X], d[X], ng);
            p2[Y] = PointSpace(j, sMin[Y], d[Y], ng);
            p2[Z] = PointSpace(k, sMin[Z], d[Z], ng);
            PrimitiveByConservative(model->gamma, model->gasR, node[idx].U[TO], Uo);
            fprintf(filePointer, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                    p2[X], p2[Y], p2[Z], Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]); 
        }
        fclose(filePointer); /* close current opened file */
    }
    return;
}
void WriteFieldDataAtCurveProbes(const Time *time, const Space *space, const Model *model)
{
    if (0 == time->curveProbeN) {
        return;
    }
    FILE *filePointer = NULL;
    String fileName = {'\0'};
    const Partition *restrict part = &(space->part);
    const Node *const node = space->node;
    const Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    int idx = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0};
    const IntVec nMin = {part->ng, part->ng, part->ng};
    const IntVec nMax = {part->n[X] - part->ng, part->n[Y] - part->ng, part->n[Z] - part->ng};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const int ng = part->ng;
    RealVec pG = {0.0}; /* ghost point */
    RealVec pO = {0.0}; /* boundary point */
    RealVec pI = {0.0}; /* image point */
    RealVec N = {0.0}; /* normal */
    int box[DIMS][LIMIT] = {{0}}; /* bounding box in node space */
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        snprintf(fileName, sizeof(fileName), "%s%03d_%05d.csv", "curve_probe_", n + 1, time->stepC);
        filePointer = fopen(fileName, "w");
        if (NULL == filePointer) {
            FatalError("failed to write data at probes...");
        }
        fprintf(filePointer, "# x, y, z, Nx, Ny, Nz, rho, u, v, w, p, T <time=%.6g>\n", time->now); 
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ValidNodeSpace(NodeSpace(poly->box[s][MIN], sMin[s], dd[s], ng), nMin[s], nMax[s]);
            box[s][MAX] = ValidNodeSpace(NodeSpace(poly->box[s][MAX], sMin[s], dd[s], ng), nMin[s], nMax[s]) + 1;
        }
        for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
            for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    if ((1 != node[idx].gst) || (n + 1 != node[idx].gid)) {
                        continue;
                    }
                    pG[X] = PointSpace(i, sMin[X], d[X], ng);
                    pG[Y] = PointSpace(j, sMin[Y], d[Y], ng);
                    pG[Z] = PointSpace(k, sMin[Z], d[Z], ng);
                    ComputeGeometricData(node[idx].fid, poly, pG, pO, pI, N);
                    PrimitiveByConservative(model->gamma, model->gasR, node[idx].U[TO], Uo);
                    fprintf(filePointer, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                            pO[X], pO[Y], pO[Z], N[X], N[Y], N[Z], Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]); 
                }
            }
        }
        fclose(filePointer); /* close current opened file */
    }
    return;
}
void WriteSurfaceForceData(const Time *time, const Space *space)
{
    if (0 == time->forceProbeN) {
        return;
    }
    FILE *filePointer = NULL;
    String fileName = {'\0'};
    const Geometry *geo = &(space->geo);
    Polyhedron *poly = NULL;
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        snprintf(fileName, sizeof(fileName), "%s%03d.csv", "surface_force_", n + 1);
        filePointer = fopen(fileName, "a");
        if (NULL == filePointer) {
            FatalError("failed to write data at probes...");
        }
        if (0 == time->stepC) { /* this is the initialization step */
            fprintf(filePointer, "# time, Fpx, Fpy, Fpz, Fvx, Fvy, Fvz, Ttx, Tty, Ttz\n"); 
        }
        fprintf(filePointer, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                time->now, poly->Fp[X], poly->Fp[Y], poly->Fp[Z], 
                poly->Fv[X], poly->Fv[Y], poly->Fv[Z], 
                poly->Tt[X], poly->Tt[Y], poly->Tt[Z]); 
        fclose(filePointer); /* close current opened file */
    }
    return;
}
/* a good practice: end file with a newline */

