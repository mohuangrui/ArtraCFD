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
    Real Uo[DIMUo] = {0.0};
    for (int n = 0; n < flow->tallyProbe; ++n) {
        snprintf(fileName, sizeof(fileName), "%s%d.%05d", "probe", n + 1, stepCount);
        filePointer = fopen(fileName, "w");
        if (NULL == filePointer) {
            FatalError("failed to write data at probes...");
        }
        fprintf(filePointer, "# points      rho     u       v       w       p\n"); 
        /* compute and adjust index range into flow region */
        int iA = FlowRegionI(ComputeI(flow->probe[n][0], space), part);
        int jA = FlowRegionJ(ComputeJ(flow->probe[n][1], space), part);
        int kA = FlowRegionK(ComputeK(flow->probe[n][2], space), part);
        int iB = FlowRegionI(ComputeI(flow->probe[n][3], space), part);
        int jB = FlowRegionJ(ComputeJ(flow->probe[n][4], space), part);
        int kB = FlowRegionK(ComputeK(flow->probe[n][5], space), part);
        int stepN = flow->probe[n][ENTRYPROBE-1] - 1;
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
            PrimitiveByConservative(Uo, idx, U, flow);
            fprintf(filePointer, "%d     %.6g      %.6g     %.6g      %.6g      %.6g    %.6g\n",
                    m, Uo[0], Uo[1], Uo[2], Uo[3], Uo[4], Uo[5]); 
        }
        fclose(filePointer); /* close current opened file */
    }
    return 0;
}
/* a good practice: end file with a newline */

