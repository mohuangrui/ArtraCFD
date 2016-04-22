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
 * Calculations are node based, boundaries are aligned with
 * the nodes. For m inner cells, there are m + 1 node layers.
 * In this program, 2D and 3D space are unified, a 2D space
 * is equivalent to a non-zero thickness 3D space with 2 inner
 * cells (that is, three node layers) in the collapsed direction.
 * These three node layers are treated as domain boundary,
 * inner node, domain boundary respectively. Zero gradient
 * condition need to be forced on the collapsed dimension.
 */
static int NodeBasedMeshNumberRefine(Space *space, const Model *model)
{
    /* set ghost layers according to numerical scheme */
    if (TVD == model->scheme) {
        space->part.ng = 1;
    }
    if (WENO == model->scheme) {
        space->part.ng = 2;
    }
    /* check and mark collapsed space. */
    space->part.collapsed = COLLAPSEN;
    if (0 == (space->part.m[Z] - 1)) {
        space->part.collapsed = COLLAPSEZ;
    }
    if (0 == (space->part.m[Y] - 1)) {
        space->part.collapsed = 2 * space->part.collapsed + COLLAPSEY;
    }
    if (0 == (space->part.m[X] - 1)) {
        space->part.collapsed = 2 * space->part.collapsed + COLLAPSEX;
    }
    for (int s = 0; s < DIMS; ++s) {
        /* ensure at least two inner cells per dimension */
        space->part.m[s] = MinInt(space->part.m[s], 2);
        /* total number of nodes need to add ghosts nodes */
        space->part.n[s] = space->part.m[s] + 1 + 2 * space->part.ng; 
    }
    space->part.totalN = space->part.n[Z] * space->part.n[Y] * space->part.n[X];
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
        space->part.d[s] = ((space->part.domain[s][MAX] - space->part.domain[s][MIN]) /
                (Real)(space->part.m[s])) / model->refL;
        space->part.domain[s][MAX] = space->part.domain[s][MAX] / model->refL;
        space->part.domain[s][MIN] = space->part.domain[s][MIN] / model->refL;
        space->part.dd[s] = 1.0 / space->part.d[s];
    }
    space->part.tinyL = 1.0e-6 * MinReal(space->part.d[X], MinReal(space->part.d[Y], space->part.d[Z]));
    /* time */
    time->end = time->end * model->refV / model->refL;
    if (0 > time->stepN) {
        time->stepN = INT_MAX;
    }
    /* model */
    if (0 > model->layers) {
        model->layers = INT_MAX;
    }
    model->gamma = 1.4;
    model->gasR = 287.058;
    /* reference Mach number */
    model->refMa = model->refV / sqrt(model->gamma * model->gasR * model->refT);
    /* reference dynamic viscosity for viscosity normalization */
    model->refMu = model->refMu / (model->refRho * model->refV * model->refL);
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

