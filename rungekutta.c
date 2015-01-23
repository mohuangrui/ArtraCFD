/****************************************************************************
 * Numeric Scheme in Time Doamin                                            *
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
#include "timer.h"
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int RungeKuttaTimeMarching(Field *field, Flux *flux, Space *space, 
        Particle *particle, Time *time, const Partition *part, const Flow *flow)
{
    ShowInformation("Time marching...");
    Real exportTimeInterval = (time->totalTime - time->currentTime);
    /* check whether current time is equal to or larger than the total time */
    if (exportTimeInterval <= 0) {
        ShowInformation("  current time is equal to or larger than total time...");
        ShowInformation("Session End");
        return 1;
    }
    /* obtain the desired export time interval */
    exportTimeInterval = exportTimeInterval / time->totalOutputTimes;
    Real accumulatedTime = 0; /* used for control when to export data */
    /* set some timers for monitoring timeconsuming of process */
    Timer operationTimer; /* timer for computing operations */
    Real operationTime = 0; /* record consuming time of operation */
    /* time marching */
    for (time->stepCount += 1; time->currentTime < time->totalTime; ++time->stepCount) {
        fprintf(stdout, "\nStep=%d; Time=%.6g\n", time->stepCount, time->currentTime);
        /*
         * Calculate dt for current time step
         */
        time->dt = exportTimeInterval;
        /*
         * Update current time stamp, if current time exceeds the total time, 
         * recompute the value of dt to make current time equal total time.
         */
        time->currentTime = time->currentTime + time->dt;
        if (time->currentTime > time->totalTime) { /* need to refine "dt" to reach totTime  */
            time->dt = time->totalTime - (time->currentTime - time->dt);
            time->currentTime = time->totalTime;
        }
        /*
         * Compute field data in current time step
         */
        /*
         * Export computed data. Use accumulatedTime as a flag, if
         * accumulatedTime increases to anticipated export interval,
         * write data out. Because the accumulatedTime very likely
         * can not increase to the exporting interval at the last
         * phase, then add a extra condition that if current time
         * is the total time, then also write data out.
         */
        accumulatedTime = accumulatedTime + time->dt;
        if ((accumulatedTime >=  exportTimeInterval) || 
                (fabs(time->currentTime - time->totalTime) < 1e-100)) {
            ++time->outputCount; /* export count increase */
            TickTime(&operationTimer);
            WriteComputedDataEnsight(field->Uo, space, particle, time, part);
            operationTime = TockTime(&operationTimer);
            accumulatedTime = 0; /* reset accumulated time */
            fprintf(stdout, "  data export time consuming: %.6lgs\n", operationTime);
        }
        /* fluid solid coupling */
        /* particle dynamics */
        /* recompute domain geometry and remeshing */
    }
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

