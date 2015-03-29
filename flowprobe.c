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
 * Static Function Declarations
 ****************************************************************************/
static int Min(const int x, const int y);
static int Max(const int x, const int y);
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
            FatalError("failed to write data at probes...");
        }
        fprintf(filePointer, "# points      rho     u       v       w       p       T\n"); 
        int iA = (int)((flow->probePos[n][0] - space->xMin) * space->ddx) + space->ng;
        int jA = (int)((flow->probePos[n][1] - space->yMin) * space->ddy) + space->ng;
        int kA = (int)((flow->probePos[n][2] - space->zMin) * space->ddz) + space->ng;
        int iB = (int)((flow->probePos[n][3] - space->xMin) * space->ddx) + space->ng;
        int jB = (int)((flow->probePos[n][4] - space->yMin) * space->ddy) + space->ng;
        int kB = (int)((flow->probePos[n][5] - space->zMin) * space->ddz) + space->ng;
        /* adjust index range into flow region */
        iA = Min(part->iSup[0] - 1, Max(part->iSub[0], iA));
        jA = Min(part->jSup[0] - 1, Max(part->jSub[0], jA));
        kA = Min(part->kSup[0] - 1, Max(part->kSub[0], kA));
        iB = Min(part->iSup[0] - 1, Max(part->iSub[0], iB));
        jB = Min(part->jSup[0] - 1, Max(part->jSub[0], jB));
        kB = Min(part->kSup[0] - 1, Max(part->kSub[0], kB));
        int stepN = flow->probe[n] - 1;
        if (1 > stepN) { /* set to lowest resolution if happens */
            stepN = 1;
        }
        if ((abs(iB - iA) < stepN) && (abs(jB - jA) < stepN) && (abs(kB - kA) < stepN)) {
            /* set to highest resolution allowed */
            stepN = Max(abs(iB - iA), Max(abs(jB - jA), abs(kB - kA)));
        }
        const Real xStep = (Real)(iB - iA) / (Real)(stepN);
        const Real yStep = (Real)(jB - jA) / (Real)(stepN);
        const Real zStep = (Real)(kB - kA) / (Real)(stepN);
        for (int m = 0; m <= stepN; ++m) {
            const int k = kA + (int)(m * zStep);
            const int j = jA + (int)(m * yStep);
            const int i = iA + (int)(m * xStep);
            idx = IndexMath(k, j, i, space) * space->dimU;
            PrimitiveByConservative(Uo, idx, U, flow);
            fprintf(filePointer, "%d     %.6g      %.6g     %.6g      %.6g      %.6g     %.6g\n",
                    m, Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]); 
        }
        fclose(filePointer); /* close current opened file */
    }
    return 0;
}
static int Min(const int x, const int y)
{
    if (x < y) {
        return x;
    }
    return y;
}
static int Max(const int x, const int y)
{
    if (x > y) {
        return x;
    }
    return y;
}
/* a good practice: end file with a newline */

