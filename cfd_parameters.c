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
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NodeBasedMeshNumberRefine(Space *, const Model *);
static int InitializeCFDParameters(Space *, Time *, Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ComputeCFDParameters(Space *space, Time *time, Model *model)
{
    ShowInformation("Computing parameters...");
    NodeBasedMeshNumberRefine(space, model);
    InitializeCFDParameters(space, time, model);
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
static int NodeBasedMeshNumberRefine(Space *space, const Model *model)
{
    /* set ghost layers according to numerical scheme */
    if (0 == model->scheme) { /* TVD */
        space->ng = 1;
    } else { /* WENO */
        space->ng = 2;
    }
    /* 
     * Check and mark collapsed space.
     * 0 -- no collapse
     * 1 -- x collapse
     * 2 -- y collapse
     * 3 -- z collapse
     * 5 -- y and x collapse
     * 7 -- z and x collapse
     * 8 -- z and y collapse
     * 17 -- z, y, and x collaspe
     */
    space->collapsed = 0; /* set to no collapse */
    if (0 == (space->nz - 1)) {
        space->collapsed = 3;
    }
    if (0 == (space->ny - 1)) {
        space->collapsed = 2 * space->collapsed + 2;
    }
    if (0 == (space->nx - 1)) {
        space->collapsed = 2 * space->collapsed + 1;
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
static int InitializeCFDParameters(Space *space, Time *time, Model *model)
{
    /* space */
    space->dz = ((space->zMax - space->zMin) / (Real)(space->nz - 1)) / model->refLength;
    space->dy = ((space->yMax - space->yMin) / (Real)(space->ny - 1)) / model->refLength;
    space->dx = ((space->xMax - space->xMin) / (Real)(space->nx - 1)) / model->refLength;
    space->zMax = space->zMax / model->refLength;
    space->yMax = space->yMax / model->refLength;
    space->xMax = space->xMax / model->refLength;
    space->zMin = space->zMin / model->refLength;
    space->yMin = space->yMin / model->refLength;
    space->xMin = space->xMin / model->refLength;
    space->ddz = 1.0 / space->dz;
    space->ddy = 1.0 / space->dy;
    space->ddx = 1.0 / space->dx;
    space->tinyL = 1.0e-6 * MinReal(space->dx, MinReal(space->dz, space->dy));
    /* time */
    time->end = time->end * model->refVelocity / model->refLength;
    if ((0 > time->stepN)) {
        time->stepN = INT_MAX;
    }
    /* fluid and flow */
    if ((0 > model->layers)) {
        model->layers = INT_MAX;
    }
    model->gamma = 1.4;
    model->gasR = 8.314462175;
    model->pi = acos(-1);
    /* reference Mach number */
    model->refMa = model->refVelocity / sqrt(model->gamma * model->gasR * model->refTemperature);
    /* reference dynamic viscosity for viscosity normalization and modify Sutherland's law */
    model->refMu = model->refMu / (model->refDensity * model->refVelocity * model->refLength);
    /*
     * Now replace some parameters by general forms that are valid
     * for both dimensional and nondimensional N-S equations, since
     * dimensional forms can be seen as normalized by reference 1.
     */
    model->gasR = 1.0 / (model->gamma * model->refMa * model->refMa);
    model->cv = model->gasR / (model->gamma - 1.0);
    return 0;
}
/* a good practice: end file with a newline */

