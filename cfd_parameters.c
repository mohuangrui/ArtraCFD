/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>                        *
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
    if (TVD == model->scheme) {
        space->ng = 1;
    }
    if (WENO == model->scheme) {
        space->ng = 2;
    }
    /* check and mark collapsed space. */
    space->collapsed = COLLAPSEN;
    if (0 == (space->m[Z] - 1)) {
        space->collapsed = COLLAPSEZ;
    }
    if (0 == (space->m[Y] - 1)) {
        space->collapsed = 2 * space->collapsed + COLLAPSEY;
    }
    if (0 == (space->m[X] - 1)) {
        space->collapsed = 2 * space->collapsed + COLLAPSEX;
    }
    for (int s = 0; s < DIMS; ++s) {
        /* change from number of cells to number of node layers */
        space->m[s] = space->m[s] + 2;
        /* total number of nodes need to add ghosts nodes */
        space->n[s] = space->m[s] + 2 * space->ng; 
    }
    space->totalN = space->n[Z] * space->n[Y] * space->n[X];
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
    for (int s = 0; s < DIMS; ++s) {
        space->d[s] = ((space->domain[s][MAX] - space->domain[s][MIN]) / (Real)(space->m[s] - 1)) / model->refLength;
        space->domain[s][MAX] = space->domain[s][MAX] / model->refLength;
        space->domain[s][MIN] = space->domain[s][MIN] / model->refLength;
        space->dd[s] = 1.0 / space->d[s];
    }
    space->tinyL = 1.0e-6 * MinReal(space->d[Z], MinReal(space->d[Y], space->d[X]));
    /* time */
    time->end = time->end * model->refVelocity / model->refLength;
    if (0 > time->stepN) {
        time->stepN = INT_MAX;
    }
    /* fluid and flow */
    if (0 > model->layers) {
        model->layers = INT_MAX;
    }
    model->gamma = 1.4;
    model->gasR = 8.314462175;
    /* reference Mach number */
    model->refMa = model->refVelocity / sqrt(model->gamma * model->gasR * model->refTemperature);
    /* reference dynamic viscosity for viscosity normalization */
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

