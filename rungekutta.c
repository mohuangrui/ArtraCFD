/****************************************************************************
 * Numeric Scheme in Time Doamin                                            *
 * Last-modified: 20 Jan 2015 11:38:57 PM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines the numeric schemes of time domain.                  *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "rungekutta.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "gcibm.h"
#include "ensight.h"
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int RungeKuttaTimeMarching(Field *field, Flux *flux, Space *space, 
        Particle *particle, Time *time, const Partition *part, const Flow *flow)
{
    ShowInformation("Time marching...");
    double exportTimeInterval = (time->totalTime - time->currentTime) / time->totalOutputTimes;
    double accumulatedTime = 0;
    /* time marching */
    for (time->stepCount += 1; time->currentTime < time->totalTime; ++time->stepCount) {
        /* calculate the dt of the current time step */
        time->dt = exportTimeInterval;
        /* update current time */
        time->currentTime = time->currentTime + time->dt;
        if (time->currentTime > time->totalTime) { /* need to refine "dt" to reach totTime  */
            time->dt = time->totalTime - (time->currentTime - time->dt);
            time->currentTime = time->totalTime;
        }
        /* compute field data */
        printf("\nStep=%d; Time=%lg; dt=%lg\n", time->stepCount, time->currentTime, time->dt);
        /* export computed data */
        accumulatedTime = accumulatedTime + time->dt;
        if ((accumulatedTime >=  exportTimeInterval) || 
                (fabs(time->currentTime - time->totalTime) < 1e-15)) {
            ++time->outputCount; /* export count increase */
            WriteComputedDataEnsight(field->Uo, space, particle, time, part);
            accumulatedTime = 0; /* reset accumulated time */
        }
        /* fluid solid coupling */
        /* particle dynamics */
        /* recompute domain geometry and remeshing */
    }
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

