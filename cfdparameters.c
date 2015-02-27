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
 * we simply use available node layers in the original cells 
 * configurations. Therefore, to have the same number of nodes layers as
 * in a cell based approach, we need to explicitly use one more cell(mesh).
 *
 * Moreover, if we have m cells without adding nodes at cell center but using
 * original nodes layers as computational domain, the total number of node
 * layers will be (m+1). To use n to stand for the total number of normal nodes
 * (that is, including boundary nodes and interior nodes but excluding exterior
 * ghost nodes), an extra 1 need to be added. As a conclusion, we need to
 * refine the inputed mesh number by 2.
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
    space->nz = space->nz + 2;
    space->ny = space->ny + 2;
    space->nx = space->nx + 2;
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
    space->dx = (space->dx / (space->nx - 1)) / flow->refLength; /* normalized real dx */
    space->dy = (space->dy / (space->ny - 1)) / flow->refLength; /* normalized real dy */
    space->dz = (space->dz / (space->nz - 1)) / flow->refLength; /* normalized real dz */
    if (space->dx <= 0) { /* zero value should have zero reciprocal */
        space->dx = 1e38;
    }
    if (space->dy <= 0) {
        space->dy = 1e38;
    }
    if (space->dz <= 0) {
        space->dz = 1e38;
    }
    /* time */
    time->totalTime = time->totalTime * flow->refVelocity / flow->refLength;
    /* fluid and flow */
    flow->gamma = 1.4;
    flow->gasR = 8.314462175;
    /* reference Mach number */
    flow->refMa = flow->refVelocity / sqrt(flow->gamma * flow->gasR * flow->refTemperature);
    /* user defined reference dynamic viscosity to simply Sutherland's law */
    flow->refMu = 1.45e-6 / (flow->refDensity * flow->refVelocity * flow->refLength);
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

