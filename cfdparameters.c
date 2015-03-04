/****************************************************************************
 * Compute and Define CFD Numerical Parameters                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file computes the numeric parameters of CFD.                      *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "cfdparameters.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NodeBasedMeshNumberRefine(Space *);
static int InitializeCFDParameters(Space *, Time *, Flow *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ComputeCFDParameters(Space *space, Time *time, Flow *flow)
{
    ShowInformation("Computing CFD parameters...");
    NodeBasedMeshNumberRefine(space);
    InitializeCFDParameters(space, time, flow);
    ShowInformation("Session End");
    return 0;
}
/*
 * Calculations in this program is node based. Instead of implicitly 
 * adding a calculation node at each cell center as in cell based 
 * program does and resulting grid size inconsistency at the boundaries,
 * we simply use available node layers in the original cell configurations. 
 * Therefore, at least two cells are required in each direction.
 *
 * Moreover, if we have m cells without adding nodes at cell center but using
 * original nodes layers as computational domain, the total number of node
 * layers will be (m+1).
 *
 * After these operations, the member n(x,y,z) in Space is the total number of
 * normal node layers. Considering the exterior ghost cells, The detailed 
 * node range of the domain will be the following:
 *
 * (In this program, Sub and Sup are used for domain range identifying,
 * therefore they are widely used in for loop control. To reduce the for loop
 * operations, we always use a reachable value for Sub, and unreachable value
 * for the Sup. To count the total valid objects, simply do Sup - Sub).
 *
 * Entire domain (includes exterior ghost): Sub = 0; Sup = n + 2 * ng;
 * Lower exterior ghost: Sub = 0; Sup = ng;
 * Normal nodes: Sub = ng; Sup = n + ng;
 *      Lower boundary: Sub = ng; Sup = ng + 1;
 *      Interior cells(node layers): Sub = ng + 1; Sup = n + ng - 1;
 *      Upper Boundary: Sub = n + ng - 1; Sup = n + ng;
 * Upper exterior ghost: Sub = n + ng; Sup = n + 2 *ng;
 *
 * In this program, 2D and 3D space are unified, that is, a 2D space
 * will automatically transfer to a zero-thickness 3D space with
 * 2 cells(that is, three node layers) in the collapsed direction.
 * These three node layers are treated as Domain Boundary, 
 * inner cell, Domain Boundary respectively. Periodic boundary 
 * condition will be forced on these two boundaries. Here the concept
 * that a zero-thickness 3D space with 2 cells is that the space is 
 * physically zero-thickness(it's still 2D), but numerically has
 * two cells(three node layers) in this direction.  
 */
static int NodeBasedMeshNumberRefine(Space *space)
{
    if (space->nz < 2) {
        space->nz = 2; /* at least two cells are required */
    }
    if (space->ny < 2) {
        space->ny = 2; /* at least two cells are required */
    }
    if (space->nx < 2) {
        space->nx = 2; /* at least two cells are required */
    }
    /* change from number of cells to number of node layers */
    space->nz = space->nz + 1;
    space->ny = space->ny + 1;
    space->nx = space->nx + 1;
    space->kMax = space->nz + 2 * space->ng; /* nz nodes + 2*ng ghosts */
    space->jMax = space->ny + 2 * space->ng; /* ny nodes + 2*ng ghosts */
    space->iMax = space->nx + 2 * space->ng; /* nx nodes + 2*ng ghosts */
    space->nMax = space->kMax * space->jMax * space->iMax; /* total node number */
    return 0;
}
/*
 * This function computes values of necessary parameters, which initializes
 * the numerical simulation environment. These parameters will be normalized
 * by reference values.
 */
static int InitializeCFDParameters(Space *space, Time *time, Flow *flow)
{
    /* space */
    space->dx = ((space->xMax - space->xMin) / (space->nx - 1)) / flow->refLength;
    space->dy = ((space->yMax - space->yMin) / (space->ny - 1)) / flow->refLength;
    space->dz = ((space->zMax - space->zMin) / (space->nz - 1)) / flow->refLength;
    space->xMax = space->xMax / flow->refLength;
    space->yMax = space->yMax / flow->refLength;
    space->zMax = space->zMax / flow->refLength;
    space->xMin = space->xMin / flow->refLength;
    space->yMin = space->yMin / flow->refLength;
    space->zMin = space->zMin / flow->refLength;
    if (space->dx > 0) {
        space->ddx = 1 / space->dx;
    } else { /* zero mesh size has zero reciprocal */
        space->ddx = 0;
    }
    if (space->dy > 0) {
        space->ddy = 1 / space->dy;
    } else { /* zero mesh size has zero reciprocal */
        space->ddy = 0;
    }
    if (space->dz > 0) {
        space->ddz = 1 / space->dz;
    } else { /* zero mesh size has zero reciprocal */
        space->ddz = 0;
    }
    /* time */
    time->totalTime = time->totalTime * flow->refVelocity / flow->refLength;
    /* fluid and flow */
    flow->gamma = 1.4;
    flow->gasR = 8.314462175;
    /* reference Mach number */
    flow->refMa = flow->refVelocity / sqrt(flow->gamma * flow->gasR * flow->refTemperature);
    /* reference dynamic viscosity for viscosity normalization */
    flow->refMu = flow->refDensity * flow->refVelocity * flow->refLength;
    /*
     * Now replace some parameters by general forms that are valid
     * for both dimensional and nondimensional N-S equations, since
     * dimensional forms can be seen as normalized by reference 1.
     */
    flow->gasR = 1 / (flow->gamma * flow->refMa * flow->refMa);
    flow->cv = flow->gasR / (flow->gamma - 1);
    return 0;
}
/* a good practice: end file with a newline */

