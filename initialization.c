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
static int ApplySpecialInitializer(const int, const Real **, Real *U, const Space *, 
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
    /* extract initial values */
    const Real rho = part->valueBC[0][0];
    const Real u = part->valueBC[0][1];
    const Real v = part->valueBC[0][2];
    const Real w = part->valueBC[0][3];
    const Real p = part->valueBC[0][4];
    /*
     * Initialize the interior field
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* linear array index math variable */
    for (k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                U[idx+0] = rho;
                U[idx+1] = rho * u;
                U[idx+2] = rho * v;
                U[idx+3] = rho * w;
                U[idx+4] = p / (flow->gamma - 1) + 0.5 * rho * (u * u + v * v + w * w);
            }
        }
    }
    /*
     * Special initializer for specific flow regions
     */
    const Real *valueIC = part->valueIC; /* pointer to the value queue of initial values */
    for (i = 1; i <= part->typeIC[0]; ++i) {
        ApplySpecialInitializer(part->typeIC[i], &valueIC, U, space, part, flow);
    }
    /*
     * Apply boundary conditions to obtain an entire initialized flow field
     */
    BoundaryCondtion(U, space, particle, part, flow);
    return 0;
}
/*
 * The handling of special initialization for specific flow region is achieved
 * through the cooperation of two data structures:
 * The typeIC array keeps a list of the types of special initialization, with
 * the first element typeIC[0] as a tally, as well as a pointer to search and 
 * manipulate the array when reading in a type entry.
 * The valueIC array is a queue data structure which stored the information of
 * the corresponding IC type. In order to keep tracking the unused data, an
 * auxiliary pointer is needed.
 * IC types and corresponding information entries:
 * 1: plane (x, y, z, normalX, normalY, normalZ, rho, u, v, w, p)
 * 2: sphere (x, y, z, r, rho, u, v, w, p)
 * 3: box (xmin, ymin, zmin, xmax, ymax, zmax, rho, u, v, w, p)
 */
static int ApplySpecialInitializer(const int typeIC, const Real **valueICPointerPointer,
        Real *U, const Space *space, const Partition *part, const Flow *flow)
{
    const Real *valueIC = *valueICPointerPointer; /* get the current valueIC pointer */
    Real rho = 0;
    Real u = 0;
    Real v = 0;
    Real w = 0;
    Real p = 0;
    Real x = 0;
    Real y = 0;
    Real z = 0;
    Real r = 0;
    Real xh = 0;
    Real yh = 0;
    Real zh = 0;
    Real normalZ = 0;
    Real normalY = 0;
    Real normalX = 0;
    /*
     * First, acquire the information data entries
     */
    x = valueIC[0];
    y = valueIC[1];
    z = valueIC[2];
    switch (typeIC) {
        case 1: /* plane */
            normalX = valueIC[3];
            normalY = valueIC[4];
            normalZ = valueIC[5];
            rho = valueIC[6];
            u = valueIC[7];
            v = valueIC[8];
            w = valueIC[9];
            p = valueIC[10];
            *valueICPointerPointer = valueIC + 11; /* update pointer of valueIC queue */
            break;
        case 2: /* sphere */
            r = valueIC[3];
            rho = valueIC[4];
            u = valueIC[5];
            v = valueIC[6];
            w = valueIC[7];
            p = valueIC[8];
            *valueICPointerPointer = valueIC + 9; /* update pointer of valueIC queue */
            break;
        case 3: /* box */
            xh = valueIC[3];
            yh = valueIC[4];
            zh = valueIC[5];
            rho = valueIC[6];
            u = valueIC[7];
            v = valueIC[8];
            w = valueIC[9];
            p = valueIC[10];
            *valueICPointerPointer = valueIC + 11; /* update pointer of valueIC queue */
            break;
        default:
            break;
    }
    /*
     * Apply initial values for nodes that meets condition
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0;
    int flag = 0; /* control flag for whether current node in the region */
    for (k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                flag = 0; /* always initialize flag to zero */
                switch (typeIC) {
                    case 1: /* plane */
                        xh = (x - (i - space->ng) * space->dx) * normalX;
                        yh = (y - (j - space->ng) * space->dy) * normalY;
                        zh = (z - (k - space->ng) * space->dz) * normalZ;
                        if ((xh + yh + zh) <= 0) { /* on the plane or on the normal direction */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case 2: /* sphere */
                        xh = (x - (i - space->ng) * space->dx);
                        yh = (y - (j - space->ng) * space->dy);
                        zh = (z - (k - space->ng) * space->dz);
                        if ((xh * xh + yh * yh + zh * zh - r * r) <= 0) { /* on or in the sphere */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case 3: /* box */
                        normalX = ((i - space->ng) * space->dx - x) * ((i - space->ng) * space->dx - xh);
                        normalY = ((j - space->ng) * space->dy - y) * ((j - space->ng) * space->dy - yh);
                        normalZ = ((k - space->ng) * space->dz - z) * ((k - space->ng) * space->dz - zh);
                        if ((normalX <= 0) || (normalY <= 0) || (normalZ <= 0)) { /* on or in the box */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    default:
                        break;
                }
                if (flag == 1) { /* current node meets the condition */
                    U[idx+0] = rho;
                    U[idx+1] = rho * u;
                    U[idx+2] = rho * v;
                    U[idx+3] = rho * w;
                    U[idx+4] = p / (flow->gamma - 1) + 0.5 * rho * (u * u + v * v + w * w);
                }
            }
        }
    }
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

