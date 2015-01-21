/****************************************************************************
 * Case Settings Checker                                                    *
 * Last-modified: 18 Jan 2015 10:27:09 PM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a checker for the case settings                      *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "casedatachecker.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function do some parameter checking
 */
int CheckCaseSettingData(const Space *space, const Time *time, 
        const Fluid *fluid, const Reference *reference)
{
    ShowInformation("Preliminary case data checking ...");
    /* space */
    if ((space->dz < 0) || (space->dy < 0) || (space->dx < 0)) {
        FatalError("Negative length values in case settings");
    }
    if ((space->nz < 1) || (space->ny < 1) || (space->nx < 1)
            || (space->ng < 1)) {
        FatalError("Too small mesh values in case settings");
    }
    /* time */
    if ((time->restart < 0) || (time->restart > 1)|| (time->totalTime <= 0)
            || (time->numCFL <= 0) || (time->totalOutputTimes < 1)) {
        FatalError("Wrong values in time section of case settings");
    }
    /* fluid */
    if ((fluid->density <= 0) || (fluid->nu <= 0) || (fluid->alpha <= 0)) {
        FatalError("Wrong values in fluid section of case settings");
    }
    /* reference */
    if ((reference->length <= 0) || (reference->density <= 0) || 
            (reference->velocity <= 0) || reference->temperature <= 0) {
        FatalError("Wrong values in reference section of case settings");
    }
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

