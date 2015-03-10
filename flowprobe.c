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
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int WriteComputedDataAtProbes(const int stepCount, const Real *U, 
        const Space *space, const Flow *flow)
{
    FILE *filePointer = NULL;
    char fileName[25] = {'\0'};
    int idx = 0; /* linear array index math variable */
    Real rho = 0; 
    Real u = 0;
    Real v = 0;
    Real w = 0;
    Real eT = 0;
    Real p = 0;
    Real T = 0;
    for (int n = 1; n <= flow->probe[0]; ++n) {
        snprintf(fileName, sizeof(fileName), "%s%d.%05d", "probe", n, stepCount);
        filePointer = fopen(fileName, "w");
        if (NULL == filePointer) {
            FatalError("failed to write data to ensight data file: ensight.***...");
        }
        fprintf(filePointer, "# points      rho     u       v       w       p       T\n"); 
        const int iA = (flow->probePos[n][0] - space->xMin) * space->ddx + space->ng;
        const int jA = (flow->probePos[n][1] - space->yMin) * space->ddy + space->ng;
        const int kA = (flow->probePos[n][2] - space->zMin) * space->ddz + space->ng;
        const int iB = (flow->probePos[n][3] - space->xMin) * space->ddx + space->ng;
        const int jB = (flow->probePos[n][4] - space->yMin) * space->ddy + space->ng;
        const int kB = (flow->probePos[n][5] - space->zMin) * space->ddz + space->ng;
        int stepN = flow->probe[n] - 1;
        if (1 > stepN) {
            stepN = 1;
        }
        const int iStep = (iB - iA) / stepN;
        const int jStep = (jB - jA) / stepN;
        const int kStep = (kB - kA) / stepN;
        for (int m = 0; m < flow->probe[n]; ++m) {
            idx = (((kA + m * kStep) * space->jMax + (jA + m * jStep)) * space->iMax + (iA + m * iStep)) * 5;
            if ((space->nMax <= idx) || (0 > idx)) {
                continue;
            }
            rho = U[idx+0];
            u = U[idx+1] / rho;
            v = U[idx+2] / rho;
            w = U[idx+3] / rho;
            eT = U[idx+4] / rho;
            p = (flow->gamma - 1) * rho * (eT - 0.5 * (u * u + v * v + w * w));
            T = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;
            fprintf(filePointer, "%d     %.6g      %.6g     %.6g      %.6g      %.6g     %.6g\n", m, rho, u, v, w, p, T); 
        }
        fclose(filePointer); /* close current opened file */
    }
    return 0;
}
/* a good practice: end file with a newline */

