/****************************************************************************
 * Numerical Solving                                                        *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines the numerical solving procedures.                    *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "solve.h"
#include <stdio.h> /* standard library for input and output */
#include "initialization.h"
#include "rungekutta.h"
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int Solve(Field *field, Flux *flux, Space *space, Particle *particle,
        Time *time, const Partition *part, const Fluid *fluid,
        const Flow *flow)
{
    InitializeFlowField(field, space, particle, time, part, flow);
    RungeKuttaTimeMarching(field, flux, space, particle, time, part, fluid, flow);
    return 0;
}
/* a good practice: end file with a newline */

