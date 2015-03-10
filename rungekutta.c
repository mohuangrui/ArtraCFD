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
    /*
     * First, save the full value of current field, since this value is required
     * for the calculation of field data at intermediate stage as well as n+1.
     * Operation can be achieved by a single loop since all data are stored by
     * linear arrays.
     */
    for (int idx = 0; idx < (space->nMax * 5); ++idx) {
        field->Un[idx] = field->U[idx];
    }
    /*
     * Then solve the flow field for a time step, updated data will be stored
     * in the same storage space of inputed data.
     */
    SpatialDiscretizationAndComputation(field->U, time->dt, field->Uswap, space, particle, part, flow);
    /*
     * Now solve the updated field data for another time step.
     */
    SpatialDiscretizationAndComputation(field->U, time->dt, field->Uswap, space, particle, part, flow);
    /*
     * Calculate the intermediate stage based on the newest updated data and
     * the original field data which is stored at first. No new storage space
     * need to be introduced since all calculations are based on a single space
     * and do not require neighbours.
     */
    const Real coeA = 3.0 / 4.0; /* caution! float operation required! */
    const Real coeB = 1.0 / 4.0; /* caution! float operation required! */
    for (int idx = 0; idx < (space->nMax * 5); ++idx) {
        field->U[idx] = coeA * field->Un[idx] + coeB * field->U[idx];
    }
    /*
     * Now solve the updated field data based on the intermediate field.
     */
    SpatialDiscretizationAndComputation(field->U, time->dt, field->Uswap, space, particle, part, flow);
    /*
     * Calculate field data at n+1
     */
    const Real coeAA = 1.0 / 3.0; /* caution! float operation required! */
    const Real coeBB = 2.0 / 3.0; /* caution! float operation required! */
    for (int idx = 0; idx < (space->nMax * 5); ++idx) {
        field->U[idx] = coeAA * field->Un[idx] + coeBB * field->U[idx];
    }
    return 0;
}
/* a good practice: end file with a newline */

