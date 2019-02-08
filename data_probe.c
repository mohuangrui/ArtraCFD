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
#include "computational_geometry.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void WritePointProbeData(const Time *time, const Space *space, const Model *model)
{
    if (0 == time->dataN[PROPT]) {
        return;
    }
    FILE *fp = NULL;
    String fname = {'\0'};
    const Partition *const part = &(space->part);
    const Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0};
    const IntVec nMin = {part->ns[PHY][X][MIN], part->ns[PHY][Y][MIN], part->ns[PHY][Z][MIN]};
    const IntVec nMax = {part->ns[PHY][X][MAX], part->ns[PHY][Y][MAX], part->ns[PHY][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const IntVec ng = {part->ng[X], part->ng[Y], part->ng[Z]};
    RealVec p1 = {0.0};
    int i = 0, j = 0, k = 0;
    for (int n = 0; n < time->dataN[PROPT]; ++n) {
        snprintf(fname, sizeof(fname), "%s%03d.csv", "point_probe_", n + 1);
        fp = Fopen(fname, "a");
        if (0 == time->stepC) { /* initialization step */
            fprintf(fp, "# time, rho, u, v, w, p, T\n");
        }
        p1[X] = time->pp[n][0];
        p1[Y] = time->pp[n][1];
        p1[Z] = time->pp[n][2];
        i = ConfineSpace(MapNode(p1[X], sMin[X], dd[X], ng[X]), nMin[X], nMax[X]);
        j = ConfineSpace(MapNode(p1[Y], sMin[Y], dd[Y], ng[Y]), nMin[Y], nMax[Y]);
        k = ConfineSpace(MapNode(p1[Z], sMin[Z], dd[Z], ng[Z]), nMin[Z], nMax[Z]);
        idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
        MapPrimitive(model->gamma, model->gasR, node[idx].U[TO], Uo);
        fprintf(fp, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                time->now, Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]);
        fclose(fp);
    }
    return;
}
void WriteLineProbeData(const Time *time, const Space *space, const Model *model)
{
    if (0 == time->dataN[PROLN]) {
        return;
    }
    FILE *fp = NULL;
    String fname = {'\0'};
    const Partition *const part = &(space->part);
    const Node *const node = space->node;
    int idx = 0; /* linear array index math variable */
    int idxOld = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0};
    const IntVec nMin = {part->ns[PHY][X][MIN], part->ns[PHY][Y][MIN], part->ns[PHY][Z][MIN]};
    const IntVec nMax = {part->ns[PHY][X][MAX], part->ns[PHY][Y][MAX], part->ns[PHY][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const IntVec ng = {part->ng[X], part->ng[Y], part->ng[Z]};
    RealVec p1 = {0.0};
    RealVec p2 = {0.0};
    RealVec dl = {0.0};
    int stepN = 0;
    int i = 0, j = 0, k = 0;
    for (int n = 0; n < time->dataN[PROLN]; ++n) {
        snprintf(fname, sizeof(fname), "%s%03d_%05d.csv", "line_probe_", n + 1, time->stepC);
        fp = Fopen(fname, "w");
        fprintf(fp, "# x, y, z, rho, u, v, w, p, T <time=%.6g>\n", time->now);
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
            i = ConfineSpace(MapNode(p1[X] + m * dl[X], sMin[X], dd[X], ng[X]), nMin[X], nMax[X]);
            j = ConfineSpace(MapNode(p1[Y] + m * dl[Y], sMin[Y], dd[Y], ng[Y]), nMin[Y], nMax[Y]);
            k = ConfineSpace(MapNode(p1[Z] + m * dl[Z], sMin[Z], dd[Z], ng[Z]), nMin[Z], nMax[Z]);
            idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
            if (idxOld == idx) {
                continue;
            }
            idxOld = idx; /* record */
            p2[X] = MapPoint(i, sMin[X], d[X], ng[X]);
            p2[Y] = MapPoint(j, sMin[Y], d[Y], ng[Y]);
            p2[Z] = MapPoint(k, sMin[Z], d[Z], ng[Z]);
            MapPrimitive(model->gamma, model->gasR, node[idx].U[TO], Uo);
            fprintf(fp, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                    p2[X], p2[Y], p2[Z], Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]);
        }
        fclose(fp);
    }
    return;
}
void WriteCurveProbeData(const Time *time, const Space *space, const Model *model)
{
    if (0 == time->dataN[PROCV]) {
        return;
    }
    FILE *fp = NULL;
    String fname = {'\0'};
    const Partition *const part = &(space->part);
    const Node *const node = space->node;
    const Geometry *const geo = &(space->geo);
    const Polyhedron *poly = NULL;
    int idx = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0};
    const IntVec nMin = {part->ns[PHY][X][MIN], part->ns[PHY][Y][MIN], part->ns[PHY][Z][MIN]};
    const IntVec nMax = {part->ns[PHY][X][MAX], part->ns[PHY][Y][MAX], part->ns[PHY][Z][MAX]};
    const RealVec sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVec d = {part->d[X], part->d[Y], part->d[Z]};
    const RealVec dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const IntVec ng = {part->ng[X], part->ng[Y], part->ng[Z]};
    RealVec pG = {0.0}; /* ghost point */
    RealVec pO = {0.0}; /* boundary point */
    RealVec pI = {0.0}; /* image point */
    RealVec N = {0.0}; /* normal */
    int box[DIMS][LIMIT] = {{0}}; /* bounding box in node space */
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        snprintf(fname, sizeof(fname), "%s%03d_%05d.csv", "curve_probe_", n + 1, time->stepC);
        fp = Fopen(fname, "w");
        fprintf(fp, "# x, y, z, Nx, Ny, Nz, rho, u, v, w, p, T <time=%.6g>\n", time->now);
        /* determine search range according to bounding box of polyhedron and valid node space */
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = ConfineSpace(MapNode(poly->box[s][MIN], sMin[s], dd[s], ng[s]), nMin[s], nMax[s]);
            box[s][MAX] = ConfineSpace(MapNode(poly->box[s][MAX], sMin[s], dd[s], ng[s]), nMin[s], nMax[s]) + 1;
        }
        for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
            for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
                for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                    idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    if ((1 != node[idx].gst) || (n + 1 != node[idx].did)) {
                        continue;
                    }
                    pG[X] = MapPoint(i, sMin[X], d[X], ng[X]);
                    pG[Y] = MapPoint(j, sMin[Y], d[Y], ng[Y]);
                    pG[Z] = MapPoint(k, sMin[Z], d[Z], ng[Z]);
                    ComputeGeometricData(pG, node[idx].fid, poly, pO, pI, N);
                    MapPrimitive(model->gamma, model->gasR, node[idx].U[TO], Uo);
                    fprintf(fp, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                            pO[X], pO[Y], pO[Z], N[X], N[Y], N[Z], Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]);
                }
            }
        }
        fclose(fp);
    }
    return;
}
void WriteSurfaceForceData(const Time *time, const Space *space, const Model *model)
{
    if (0 == time->dataN[PROFC]) {
        return;
    }
    FILE *fp = NULL;
    String fname = {'\0'};
    const Geometry *const geo = &(space->geo);
    const Polyhedron *poly = NULL;
    for (int n = 0; n < geo->totN; ++n) {
        poly = geo->poly + n;
        snprintf(fname, sizeof(fname), "%s%03d.csv", "surface_force_", n + 1);
        fp = Fopen(fname, "a");
        if (0 == time->stepC) { /* initialization step */
            fprintf(fp, "# time, Fpx, Fpy, Fpz, Fvx, Fvy, Fvz, Ttx, Tty, Ttz <model.mid=%d>\n", model->mid);
        }
        fprintf(fp, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                time->now, poly->Fp[X], poly->Fp[Y], poly->Fp[Z],
                poly->Fv[X], poly->Fv[Y], poly->Fv[Z],
                poly->Tt[X], poly->Tt[Y], poly->Tt[Z]);
        fclose(fp);
    }
    return;
}
/* a good practice: end file with a newline */

