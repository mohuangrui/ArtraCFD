/****************************************************************************
 * Flow Initialization                                                      *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a function that handles flow initialization.         *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "initialization.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "boundarycondition.h"
#include "ensight.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int FirstRunInitializer(Real *U, const Space *, const Particle *,
        const Partition *, const Flow *);
static int RestartInitializer(Real *U, const Space *, Time *, 
        const Partition *, const Flow *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function initializes the entire flow field. Initialization will be 
 * done differently determined by the restart status.
 */
int InitializeFlowField(Real *U, const Space *space, const Particle *particle,
        Time *time, const Partition *part, const Flow *flow)
{
    ShowInformation("Initializing flow field...");
    if (time->restart == 0) { /* non restart */
        FirstRunInitializer(U, space, particle, part, flow);
        /* if this is a first run, output initial data */
        InitializeEnsightTransientCaseFile(time);
        WriteComputedDataEnsight(U, space, particle, time, part, flow);
    } else {
        RestartInitializer(U, space, time, part, flow);
    }
    ShowInformation("Session End");
    return 0;
}
/*
 * The first run initialization will assign values to field variables.
 */
static int FirstRunInitializer(Real *U, const Space *space, const Particle *particle, 
        const Partition *part, const Flow *flow)
{
    ShowInformation("  Non-restart run initializing...");
    /*
     * Initial conditions, these values should be
     * normalized values relative to the reference values.
     */
    const Real rho0 = 1.0;
    const Real u0 = 0.0;
    const Real v0 = 0.0;
    const Real w0 = 0.0;
    const Real p0 = 1.0;
    /*
     * Initialize the interior field
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* linear array index math variable */
    for (k = part->kSub[12]; k < part->kSup[12]; ++k) {
        for (j = part->jSub[12]; j < part->jSup[12]; ++j) {
            for (i = part->iSub[12]; i < part->iSup[12]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                U[idx+0] = rho0;
                U[idx+1] = rho0 * u0;
                U[idx+2] = rho0 * v0;
                U[idx+3] = rho0 * w0;
                U[idx+4] = p0 / (flow->gamma - 1) + 
                    0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0);
            }
        }
    }
    /*
     * Apply boundary conditions to obtain an entire initialized flow field
     */
    BoundaryCondtion(U, space, particle, part, flow);
    return 0;
}
/*
 * If this is a restart run, then initialize flow field by reading field data
 * from restart files.
 */
static int RestartInitializer(Real *U, const Space *space, Time *time,
        const Partition *part, const Flow *flow)
{
    ShowInformation("  Restart run initializing...");
    /*
     * Load data from Ensight restart files.
     */
    LoadComputedDataEnsight(U, space, time, part, flow);
    return 0;
}
/* a good practice: end file with a newline */

