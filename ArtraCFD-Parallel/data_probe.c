/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) 2014-2018 Huangrui Mo <huangrui.mo@gmail.com>              *
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
int WriteComputedDataAtProbes(const Real *U, const Space *space, 
        const Time *time, const Model *model, const Partition *part)
{
    FILE *filePointer = NULL;
    char fileName[25] = {'\0'};
    int idx = 0; /* linear array index math variable */
    Real Uo[DIMUo] = {0.0};
    for (int n = 0; n < time->tallyProbe; ++n) {
        snprintf(fileName, sizeof(fileName), "%s%d.%05d", "probe", n + 1, time->stepCount);
        filePointer = fopen(fileName, "w");
        if (NULL == filePointer) {
            FatalError("failed to write data at probes...");
        }
        fprintf(filePointer, "# points      rho     u       v       w       p       T\n"); 
        /* compute and adjust index range into flow region */
        int iA = ValidRegionI(ComputeI(time->probe[n][0], space), part);
        int jA = ValidRegionJ(ComputeJ(time->probe[n][1], space), part);
        int kA = ValidRegionK(ComputeK(time->probe[n][2], space), part);
        int iB = ValidRegionI(ComputeI(time->probe[n][3], space), part);
        int jB = ValidRegionJ(ComputeJ(time->probe[n][4], space), part);
        int kB = ValidRegionK(ComputeK(time->probe[n][5], space), part);
        int stepN = time->probe[n][ENTRYPROBE-1] - 1;
        if (1 > stepN) { /* set to lowest resolution if happens */
            stepN = 1;
        }
        if ((abs(iB - iA) < stepN) && (abs(jB - jA) < stepN) && (abs(kB - kA) < stepN)) {
            /* set to highest resolution allowed */
            stepN = MaxInt(abs(iB - iA), MaxInt(abs(jB - jA), abs(kB - kA)));
        }
        const Real xStep = (Real)(iB - iA) / (Real)(stepN);
        const Real yStep = (Real)(jB - jA) / (Real)(stepN);
        const Real zStep = (Real)(kB - kA) / (Real)(stepN);
        for (int m = 0; m <= stepN; ++m) {
            const int k = kA + (int)(m * zStep);
            const int j = jA + (int)(m * yStep);
            const int i = iA + (int)(m * xStep);
            idx = IndexMath(k, j, i, space) * DIMU;
            PrimitiveByConservative(Uo, idx, U, model);
            fprintf(filePointer, "%d     %.6g      %.6g     %.6g      %.6g      %.6g      %.6g\n",
                    m, Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]); 
        }
        fclose(filePointer); /* close current opened file */
    }
    return 0;
}
/* a good practice: end file with a newline */

