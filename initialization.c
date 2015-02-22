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
#include "cfdcommons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int FirstRunInitializer(Field *, Space *, const Partition *, const Flow *);
static int RestartInitializer(Field *, Space *, Time *, const Partition *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function initializes the entire flow field. Initialization will be 
 * done differently determined by the restart status.
 */
int InitializeFlowField(Field *field, Space *space, const Particle *particle,
        Time *time, const Partition *part, const Flow *flow)
{
    ShowInformation("Initializing flow field...");
    if (time->restart == 0) { /* non restart */
        FirstRunInitializer(field, space, part, flow);
        /* if this is a first run, output initial data */
        InitializeEnsightTransientCaseFile(time);
        WriteComputedDataEnsight(field->Uo, space, particle, time, part);
    } else {
        RestartInitializer(field, space, time,  part);
    }
    ShowInformation("Session End");
    return 0;
}
/*
 * The first run initialization will assign values to field variables.
 */
static int FirstRunInitializer(Field *field, Space *space, const Partition *part, const Flow *flow)
{
    ShowInformation("  Non-restart run initializing...");
    /*
     * Initial conditions, these values should be
     * normalized values relative to the reference values.
     */
    Real rho0 = 1.0;
    Real u0 = 0.0;
    Real v0 = 0.0;
    Real w0 = 0.0;
    Real p0 = 1.0;
    /*
     * Decompose the field variable into each component.
     */
    Real *Un[5] = {
        field->Un + 0 * space->nMax,
        field->Un + 1 * space->nMax,
        field->Un + 2 * space->nMax,
        field->Un + 3 * space->nMax,
        field->Un + 4 * space->nMax};
    /*
     * Initialize the interior field
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* calculated index */
    for (k = part->kSub[12]; k < part->kSup[12]; ++k) {
        for (j = part->jSub[12]; j < part->jSup[12]; ++j) {
            for (i = part->iSub[12]; i < part->iSup[12]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                Un[0][idx] = rho0;
                Un[1][idx] = rho0 * u0;
                Un[2][idx] = rho0 * v0;
                Un[3][idx] = rho0 * w0;
                Un[4][idx] = p0 / (flow->gamma - 1) + 
                    0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0);
            }
        }
    }
    /*
     * Apply boundary conditions to obtain an entire initialized flow field
     */
    BoundaryCondtion(field, space, part, flow);
    /*
     * Compute primitive variables based on conservative variables
     */
    ComputePrimitiveByConservative(field, space, flow);
    return 0;
}
/*
 * If this is a restart run, then initialize flow field by reading field data
 * from restart files.
 */
static int RestartInitializer(Field *field, Space *space, Time *time,
        const Partition *part)
{
    ShowInformation("  Restart run initializing...");
    /*
     * Load data from Ensight restart files.
     */
    LoadComputedDataEnsight(field->Uo, space, time, part);
    return 0;
}
/* a good practice: end file with a newline */

