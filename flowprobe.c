/****************************************************************************
 * Probes                                                                   *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines the point and line probes for flow information.      *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* support for abs operation */
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteComputedDataAtProbes(const int stepCount, const Real *U, 
        const Space *space, const Partition *part, const Flow *flow)
{
    FILE *filePointer = NULL;
    char fileName[25] = {'\0'};
    int idx = 0; /* linear array index math variable */
    Real Uo[6] = {0.0};
    for (int n = 1; n <= flow->probe[0]; ++n) {
        snprintf(fileName, sizeof(fileName), "%s%d.%05d", "probe", n, stepCount);
        filePointer = fopen(fileName, "w");
        if (NULL == filePointer) {
            FatalError("failed to write data to ensight data file: ensight.***...");
        }
        fprintf(filePointer, "# points      rho     u       v       w       p       T\n"); 
        /* plus one to shift away from boundary when dimension collapse */
        const int iA = (int)((flow->probePos[n][0] - space->xMin) * space->ddx) + space->ng;
        const int jA = (int)((flow->probePos[n][1] - space->yMin) * space->ddy) + space->ng;
        const int kA = (int)((flow->probePos[n][2] - space->zMin) * space->ddz) + space->ng;
        const int iB = (int)((flow->probePos[n][3] - space->xMin) * space->ddx) + space->ng;
        const int jB = (int)((flow->probePos[n][4] - space->yMin) * space->ddy) + space->ng;
        const int kB = (int)((flow->probePos[n][5] - space->zMin) * space->ddz) + space->ng;
        int stepN = flow->probe[n] - 1;
        if (1 > stepN) { /* set to lowest resolution if happens */
            stepN = 1;
        }
        if ((abs(iB - iA) < stepN) && (abs(jB - jA) < stepN) && (abs(kB - kA) < stepN)) {
            /* set to highest resolution allowed */
            stepN = abs(iB - iA);
            if (abs(jB - jA) > stepN) {
                stepN = abs(jB - jA);
            }
            if (abs(kB - kA) > stepN) {
                stepN = abs(kB - kA);
            }
        }
        const Real xStep = (Real)(iB - iA) / (Real)(stepN);
        const Real yStep = (Real)(jB - jA) / (Real)(stepN);
        const Real zStep = (Real)(kB - kA) / (Real)(stepN);
        for (int m = 0; m <= stepN; ++m) {
            const int k = kA + (int)(m * zStep);
            const int j = jA + (int)(m * yStep);
            const int i = iA + (int)(m * xStep);
            if ((part->kSub[0] > k) || (part->kSup[0] <= k) || 
                    (part->jSub[0] > j) || (part->jSup[0] <= j) ||
                    (part->iSub[0] > i) || (part->iSup[0] <= i)) {
                continue; /* not in flow domain */
            }
            idx = IndexMath(k, j, i, space) * space->dimU;
            PrimitiveByConservative(Uo, idx, U, flow);
            fprintf(filePointer, "%d     %.6g      %.6g     %.6g      %.6g      %.6g     %.6g\n",
                    m, Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]); 
        }
        fclose(filePointer); /* close current opened file */
    }
    return 0;
}
/* a good practice: end file with a newline */

