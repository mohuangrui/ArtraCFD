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
static int InitializeCFDParameters(Time *, Space *, Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ComputeCFDParameters(Time *time, Space *space, Model *model)
{
    NodeBasedMeshNumberRefine(space, model);
    InitializeCFDParameters(time, space, model);
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
    Partition * part = &(space->part);
    /* set ghost layers according to numerical scheme */
    if (WENOFIVE == model->scheme) {
        part->ng = 2;
    }
    /* check and mark collapsed space. */
    part->collapsed = COLLAPSEN;
    if (0 == (part->m[Z] - 1)) {
        part->collapsed = COLLAPSEZ;
    }
    if (0 == (part->m[Y] - 1)) {
        part->collapsed = 2 * part->collapsed + COLLAPSEY;
    }
    if (0 == (part->m[X] - 1)) {
        part->collapsed = 2 * part->collapsed + COLLAPSEX;
    }
    for (int s = 0; s < DIMS; ++s) {
        /* ensure at least two inner cells per dimension */
        part->m[s] = MaxInt(part->m[s], 2);
        /* total number of nodes (including ghost nodes) */
        part->n[s] = part->m[s] + 1 + 2 * part->ng; 
    }
    return 0;
}
/*
 * This function computes values of necessary parameters, which initializes
 * the numerical simulation environment. These parameters will be normalized
 * by reference values.
 */
static int InitializeCFDParameters(Time *time, Space *space, Model *model)
{
    Partition * part = &(space->part);
    /* space */
    for (int s = 0; s < DIMS; ++s) {
        part->d[s] = ((part->domain[s][MAX] - part->domain[s][MIN]) /
                (Real)(part->m[s])) / model->refL;
        part->domain[s][MAX] = part->domain[s][MAX] / model->refL;
        part->domain[s][MIN] = part->domain[s][MIN] / model->refL;
        part->dd[s] = 1.0 / part->d[s];
    }
    part->tinyL = 1.0e-6 * MinReal(part->d[X], MinReal(part->d[Y], part->d[Z]));
    /* time */
    time->end = time->end * model->refV / model->refL;
    if (0 > time->stepN) {
        time->stepN = INT_MAX;
    }
    time->countOutput = time->restart;
    /* model */
    if (0 > model->layers) {
        model->layers = INT_MAX;
    }
    model->gamma = 1.4;
    model->gasR = 287.058;
    for (int s = 0; s < DIMS; ++s) {
        model->g[s] = model->g[s] * model->refL / (model->refV * model->refV);
    }
    model->sState = model->gState; /* source state on if gravity on */
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

