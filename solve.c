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
#include "temperalmarching.h"
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int Solve(Field *field, Space *space, Particle *particle, Time *time, 
        const Partition *part, const Flow *flow)
{
    InitializeFlowField(field->Un, space, particle, time, part, flow);
    TemperalMarching(field, space, particle, time, part, flow);
    return 0;
}
/* a good practice: end file with a newline */

