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
#include <math.h> /* common mathematical functions */
#include "gcibm.h"
#include "ensight.h"
#include "timer.h"
#include "tvd.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static Real ComputeTimeStepByCFL(const Real *U, const Space *, const Time *, 
        const Partition *, const Flow *);
static Real MinPositive(const Real valueA, const Real valueB);
static Real Max(const Real valueA, const Real valueB);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int RungeKuttaTimeMarching(Field *field, Space *space, Particle *particle,
        Time *time, const Partition *part, const Flow *flow)
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
    /* set some timers for monitoring time consuming of process */
    Timer operationTimer; /* timer for computing operations */
    Real operationTime = 0; /* record consuming time of operation */
    /* time marching */
    for (time->stepCount += 1; (time->currentTime < time->totalTime) && 
            (time->stepCount <= time->stepCount); ++(time->stepCount)) {
        /*
         * Calculate dt for current time step
         */
        time->dt = ComputeTimeStepByCFL(field->Un, space, time, part, flow);
        fprintf(stdout, "\nStep=%d; Time=%.6g; dt=%.6g\n", time->stepCount, time->currentTime, time->dt);
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
                (fabs(time->currentTime - time->totalTime) < 1e-38)) {
            ++(time->outputCount); /* export count increase */
            TickTime(&operationTimer);
            WriteComputedDataEnsight(field->Un, space, particle, time, part, flow);
            operationTime = TockTime(&operationTimer);
            accumulatedTime = 0; /* reset accumulated time */
            fprintf(stdout, "  data export time consuming: %.6gs\n", operationTime);
        }
        /* fluid solid coupling */
        /* particle dynamics */
        /* recompute domain geometry and remeshing */
    }
    ShowInformation("Session End");
    return 0;
}
static Real ComputeTimeStepByCFL(const Real *U, const Space *space, const Time *time, 
        const Partition *part, const Flow *flow)
{
    /*
     * Define the primitive field variables.
     */
    Real rho = 0; 
    Real u = 0;
    Real v = 0;
    Real w = 0;
    Real p = 0;
    Real eT = 0;
    /*
     * Auxiliary variables
     */
    Real velocity = 0;
    Real velocityMax = 1e-38;
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* linear array index math variable */
    for (k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (space->ghostFlag[idx] == -1) { /* if it's solid node */
                    continue;
                }
                idx = idx * 5; /* change idx to field variable */
                rho = U[idx+0];
                u = U[idx+1] / rho;
                v = U[idx+2] / rho;
                w = U[idx+3] / rho;
                eT = U[idx+4] / rho;
                p = (flow->gamma - 1) * rho * (eT - 0.5 * (u * u + v * v + w * w));

                velocity = sqrt(flow->gamma * p / rho) + Max(fabs(u), Max(fabs(v), fabs(w)));
                if (velocityMax < velocity) {
                    velocityMax = velocity;
                }
            }
        }
    }
    return time->numCFL * MinPositive(space->dx, MinPositive(space->dy, space->dz)) / velocityMax;
}
static Real MinPositive(const Real valueA, const Real valueB)
{
    if ((valueA <= 0) && (valueB <= 0)) {
        return 1e10;
    }
    if (valueA <= 0) {
        return valueB;
    }
    if (valueB <= 0) {
        return valueA;
    }
    if (valueA < valueB) {
        return valueA;
    }
    return valueB;
}
static Real Max(const Real valueA, const Real valueB)
{
    if (valueA > valueB) {
        return valueA;
    }
    return valueB;
}
/* a good practice: end file with a newline */

