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
    Real p0 = 0.71429;
    /*
     * Decompose the field variable into each component.
     */
    Real *rho = field->Un + 0 * space->kMax * space->jMax * space->iMax;
    Real *rho_u = field->Un + 1 * space->kMax * space->jMax * space->iMax;
    Real *rho_v = field->Un + 2 * space->kMax * space->jMax * space->iMax;
    Real *rho_w = field->Un + 3 * space->kMax * space->jMax * space->iMax;
    Real *rho_eT = field->Un + 4 * space->kMax * space->jMax * space->iMax;
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
                rho[idx] = rho0;
                rho_u[idx] = rho0 * u0;
                rho_v[idx] = rho0 * v0;
                rho_w[idx] = rho0 * w0;
                rho_eT[idx] = p0 / (flow->gamma - 1) + 
                    0.5 * rho0 * (u0 * u0 + v0 * v0 + w0 * w0);
            }
        }
    }
    /*
     * Apply boundary conditions to obtain an entire initialized flow field
     */
    BoundaryCondtion(field, space, part, flow);
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

