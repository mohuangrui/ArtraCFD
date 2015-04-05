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
#include "datastream.h"
#include "geometrystream.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int FirstRunInitializer(Real *U, const Space *, const Particle *,
        const Partition *, const Flow *);
static int ApplyRegionalInitializer(const int, Real *U, const Space *, 
        const Partition *, const Flow *);
static int RestartInitializer(Real *U, const Space *, const Particle *, Time *, 
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
    if (0 == time->restart) { /* non restart */
        FirstRunInitializer(U, space, particle, part, flow);
        /* if this is a first run, output initial data */
        WriteComputedData(U, space, time, part, flow);
        WriteGeometryData(particle, time);
    } else {
        RestartInitializer(U, space, particle, time, part, flow);
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
    const Real Uo[DIMUo] = {
        part->valueBC[0][0],
        part->valueBC[0][1],
        part->valueBC[0][2],
        part->valueBC[0][3],
        part->valueBC[0][4]};
    /*
     * Initialize the interior field
     */
    int idx = 0; /* linear array index math variable */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space) * DIMU;
                ConservativeByPrimitive(U, idx, Uo, flow);
            }
        }
    }
    /*
     * Regional initializer for specific flow regions
     */
    for (int n = 0; n < part->tallyIC; ++n) {
        ApplyRegionalInitializer(n, U, space, part, flow);
    }
    /*
     * Boundary conditions and treatments to obtain an entire initialized flow field
     */
    BoundaryCondtionsAndTreatments(U, space, particle, part, flow);
    return 0;
}
/*
 * The handling of regional initialization for specific flow region is achieved
 * through the cooperation of three data structures:
 * The tallyIC counts the number of initializers.
 * The typeIC array keeps a list of the types of regional initialization.
 * The valueIC array stored the information of the corresponding IC type.
 * IC types and corresponding information entries:
 * 1: plane (x, y, z, normalX, normalY, normalZ, rho, u, v, w, p)
 * 2: sphere (x, y, z, r, rho, u, v, w, p)
 * 3: box (xmin, ymin, zmin, xmax, ymax, zmax, rho, u, v, w, p)
 */
static int ApplyRegionalInitializer(const int n, Real *U, const Space *space, 
        const Partition *part, const Flow *flow)
{
    /*
     * Acquire the specialized information data entries
     */
    /* the fix index part */
    const Real x = part->valueIC[n][0];
    const Real y = part->valueIC[n][1];
    const Real z = part->valueIC[n][2];
    const Real Uo[DIMUo] = {
        part->valueIC[n][ENTRYIC-5],
        part->valueIC[n][ENTRYIC-4],
        part->valueIC[n][ENTRYIC-3],
        part->valueIC[n][ENTRYIC-2],
        part->valueIC[n][ENTRYIC-1]};
    /* the vary part */
    Real r = 0.0;
    Real xh = 0.0;
    Real yh = 0.0;
    Real zh = 0.0;
    Real normalZ = 0.0;
    Real normalY = 0.0;
    Real normalX = 0.0;
    switch (part->typeIC[n]) {
        case 1: /* plane */
            normalX = part->valueIC[n][3];
            normalY = part->valueIC[n][4];
            normalZ = part->valueIC[n][5];
            break;
        case 2: /* sphere */
            r = part->valueIC[n][3];
            break;
        case 3: /* box */
            xh = part->valueIC[n][3];
            yh = part->valueIC[n][4];
            zh = part->valueIC[n][5];
            break;
        default:
            break;
    }
    /*
     * Apply initial values for nodes that meets condition
     */
    int idx = 0;
    int flag = 0; /* control flag for whether current node in the region */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                flag = 0; /* always initialize flag to zero */
                switch (part->typeIC[n]) {
                    case 1: /* plane */
                        xh = (ComputeX(i, space) - x) * normalX;
                        yh = (ComputeY(j, space) - y) * normalY;
                        zh = (ComputeZ(k, space) - z) * normalZ;
                        if (0 <= (xh + yh + zh)) { /* on the normal direction or the plane */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case 2: /* sphere */
                        xh = (ComputeX(i, space) - x);
                        yh = (ComputeY(j, space) - y);
                        zh = (ComputeZ(k, space) - z);
                        if (0 >= (xh * xh + yh * yh + zh * zh - r * r)) { /* in or on the sphere */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case 3: /* box */
                        xh = (ComputeX(i, space) - x) * (ComputeX(i, space) - xh);
                        yh = (ComputeY(j, space) - y) * (ComputeY(j, space) - yh);
                        zh = (ComputeZ(k, space) - z) * (ComputeZ(k, space) - zh);
                        if ((0 >= xh) && (0 >= yh) && (0 >= zh)) { /* in or on the box */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    default:
                        break;
                }
                if (1 == flag) { /* current node meets the condition */
                    idx = IndexMath(k, j, i, space) * DIMU;
                    ConservativeByPrimitive(U, idx, Uo, flow);
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
static int RestartInitializer(Real *U, const Space *space, const Particle *particle, 
        Time *time, const Partition *part, const Flow *flow)
{
    ShowInformation("  Restart run initializing...");
    /*
     * Load data from Ensight restart files.
     */
    LoadComputedData(U, space, time, part, flow);
    /*
     * Boundary conditions and treatments to obtain an entire initialized flow field
     */
    BoundaryCondtionsAndTreatments(U, space, particle, part, flow);
    return 0;
}
/* a good practice: end file with a newline */

