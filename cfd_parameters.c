/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) 2014-2018 Huangrui Mo <huangrui.mo@gmail.com>              *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "cfd_parameters.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <limits.h> /* sizes of integral types */
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
    ShowInformation("Computing parameters...");
    NodeBasedMeshNumberRefine(space);
    InitializeCFDParameters(space, time, flow);
    ShowInformation("Session End");
    return 0;
}
/*
 * Calculations in this program is node based, we first construct the
 * distribution of nodes with the first node and last node placing at
 * boundaries. Then cells are constructed based on the node distribution.
 * Therefore, boundaries are aligned with the nodes, or say, cell centers.
 *
 * Consequently, if we have m cells, we will have m + 2 node layers.
 *
 * After these operations, member nx, ny, nz in Space is the total number of
 * node layers. Considering the exterior ghost cells, The detailed 
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
 * will be equivalent to a non-zero thickness 3D space with
 * 1 cells(that is, three node layers) in the collapsed direction.
 * These three node layers are treated as Domain Boundary, 
 * inner node, Domain Boundary respectively. Periodic boundary 
 * condition need to be forced on these two boundaries.
 * That is, dimension collapse should be achieved by single cell 
 * with periodic boundary conditions;
 */
static int NodeBasedMeshNumberRefine(Space *space)
{
    /* check whether space collapsed */
    if (0 == (space->nz - 1) * (space->ny - 1) * (space->nx - 1)) {
        space->collapsed = 1;
    }
    /* change from number of cells to number of node layers */
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
    space->dz = ((space->zMax - space->zMin) / (Real)(space->nz - 1)) / flow->refLength;
    space->dy = ((space->yMax - space->yMin) / (Real)(space->ny - 1)) / flow->refLength;
    space->dx = ((space->xMax - space->xMin) / (Real)(space->nx - 1)) / flow->refLength;
    space->zMax = space->zMax / flow->refLength;
    space->yMax = space->yMax / flow->refLength;
    space->xMax = space->xMax / flow->refLength;
    space->zMin = space->zMin / flow->refLength;
    space->yMin = space->yMin / flow->refLength;
    space->xMin = space->xMin / flow->refLength;
    space->ddz = 1.0 / space->dz;
    space->ddy = 1.0 / space->dy;
    space->ddx = 1.0 / space->dx;
    space->tinyL = 1.0e-6 * MinReal(space->dx, MinReal(space->dz, space->dy));
    /* time */
    time->totalTime = time->totalTime * flow->refVelocity / flow->refLength;
    if ((0 > time->totalStep)) {
        time->totalStep = INT_MAX;
    }
    /* fluid and flow */
    flow->gamma = 1.4;
    flow->gasR = 8.314462175;
    flow->pi = acos(-1);
    /* reference Mach number */
    flow->refMa = flow->refVelocity / sqrt(flow->gamma * flow->gasR * flow->refTemperature);
    /* reference dynamic viscosity for viscosity normalization and modify Sutherland's law */
    flow->refMu = flow->refMu / (flow->refDensity * flow->refVelocity * flow->refLength);
    /*
     * Now replace some parameters by general forms that are valid
     * for both dimensional and nondimensional N-S equations, since
     * dimensional forms can be seen as normalized by reference 1.
     */
    flow->gasR = 1.0 / (flow->gamma * flow->refMa * flow->refMa);
    flow->cv = flow->gasR / (flow->gamma - 1.0);
    return 0;
}
/* a good practice: end file with a newline */

