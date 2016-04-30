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
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteFieldDataAtProbes(const Time *time, const Space *space, const Model *model)
{
    FILE *filePointer = NULL;
    String fileName = {'\0'};
    const Partition *restrict part = &(space->part);
    const Node *node = space->node;
    const Real *restrict U = NULL;
    int idx = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0};
    const IntVector nMin = {part->ns[PIN][X][MIN], part->ns[PIN][Y][MIN], part->ns[PIN][Z][MIN]};
    const IntVector nMax = {part->ns[PIN][X][MAX], part->ns[PIN][Y][MAX], part->ns[PIN][Z][MAX]};
    const RealVector sMin = {part->domain[X][MIN], part->domain[Y][MIN], part->domain[Z][MIN]};
    const RealVector dd = {part->dd[X], part->dd[Y], part->dd[Z]};
    const int ng = part->ng;
    IntVector p1 = {0};
    IntVector p2 = {0};
    RealVector step = {0.0};
    for (int n = 0; n < time->probeN; ++n) {
        snprintf(fileName, sizeof(fileName), "%s%d.%05d", "probe", n + 1, time->countStep);
        filePointer = fopen(fileName, "w");
        if (NULL == filePointer) {
            FatalError("failed to write data at probes...");
        }
        fprintf(filePointer, "# points      rho     u       v       w       p       T\n"); 
        /* compute and adjust index range into flow region */
        p1[X] = ValidNodeSpace(NodeSpace(time->probe[n][0], sMin[X], dd[X], ng), nMin[X], nMax[X]);
        p1[Y] = ValidNodeSpace(NodeSpace(time->probe[n][1], sMin[Y], dd[Y], ng), nMin[Y], nMax[Y]);
        p1[Z] = ValidNodeSpace(NodeSpace(time->probe[n][2], sMin[Z], dd[Z], ng), nMin[Z], nMax[Z]);
        p2[X] = ValidNodeSpace(NodeSpace(time->probe[n][3], sMin[X], dd[X], ng), nMin[X], nMax[X]);
        p2[Y] = ValidNodeSpace(NodeSpace(time->probe[n][4], sMin[Y], dd[Y], ng), nMin[Y], nMax[Y]);
        p2[Z] = ValidNodeSpace(NodeSpace(time->probe[n][5], sMin[Z], dd[Z], ng), nMin[Z], nMax[Z]);
        int stepN = time->probe[n][ENTRYPROBE-1] - 1;
        if (1 > stepN) { /* set to lowest resolution if happens */
            stepN = 1;
        }
        if ((abs(p2[X] - p1[X]) < stepN) && (abs(p2[Y] - p1[Y]) < stepN) && (abs(p2[Z] - p1[Z]) < stepN)) {
            /* set to highest resolution allowed */
            stepN = MaxInt(abs(p2[X] - p1[X]), MaxInt(abs(p2[Y] - p1[Y]), abs(p2[Z] - p1[Z])));
        }
        step[X] = (Real)(p2[X] - p1[X]) / (Real)(stepN);
        step[Y] = (Real)(p2[Y] - p1[Y]) / (Real)(stepN);
        step[Z] = (Real)(p2[Z] - p1[Z]) / (Real)(stepN);
        for (int m = 0; m <= stepN; ++m) {
            const int k = p1[X] + (int)(m * step[X]);
            const int j = p1[Y] + (int)(m * step[Y]);
            const int i = p1[Z] + (int)(m * step[Z]);
            idx = IndexNode(k, j, i, part->n[Y], part->n[X]);
            U = node[idx].U[TO];
            PrimitiveByConservative(model->gamma, model->gasR, U, Uo);
            fprintf(filePointer, "%d     %.6g      %.6g     %.6g      %.6g      %.6g      %.6g\n",
                    m, Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]); 
        }
        fclose(filePointer); /* close current opened file */
    }
    return 0;
}
/* a good practice: end file with a newline */

