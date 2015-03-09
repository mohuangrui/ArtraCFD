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
#include "tvd.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int RungeKutta(Field *field, Space *space, Particle *particle,
        Time *time, const Partition *part, const Flow *flow)
{
    SpatialDiscretizationAndComputation(field->U, field->Un, space, particle, part, flow, time->dt);
    return 0;
}
/* a good practice: end file with a newline */

