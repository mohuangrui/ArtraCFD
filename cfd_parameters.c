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
static void SetNodeNumber(Space *, Model *);
static void InitializeParameters(Time *, Space *, Model *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void ComputeParameters(Time *time, Space *space, Model *model)
{
    SetNodeNumber(space, model);
    InitializeParameters(time, space, model);
    return;
}
/*
 * Calculations are node based, global domain boundaries are aligned
 * with node layers. On each dimension, for m inner cells, there are
 * m + 1 node layers. 2D and 3D space are unified, a 2D space is
 * equivalent to a non-zero thickness 3D space with two inner cells
 * (three node layers) in the collapsed direction. These three node
 * layers are treated as domain boundary, inner node, domain boundary
 * respectively. Zero gradient condition need to be enforced on the
 * collapsed dimensions. Using three rather than one node layers is to
 * be compatible with the three-dimensional governing equations,
 * especially the calculation of the diffusive fluxes.
 */
static void SetNodeNumber(Space *space, Model *model)
{
    Partition *const part = &(space->part);
    /* check and mark collapsed space */
    part->collapse = COLLAPSEN;
    if (0 == (part->m[Z] - 1)) {
        part->collapse = COLLAPSEZ;
    }
    if (0 == (part->m[Y] - 1)) {
        part->collapse = 2 * part->collapse + COLLAPSEY;
    }
    if (0 == (part->m[X] - 1)) {
        part->collapse = 2 * part->collapse + COLLAPSEX;
    }
    /* set stencil width and ghost layers required by numerical scheme */
    switch (model->sScheme) {
        case WENOTHREE:
            model->sL = -1; model->sR = 2; part->gl = 2;
            break;
        case WENOFIVE:
            model->sL = -2; model->sR = 3; part->gl = 3;
            break;
        default:
            break;
    }
    /*
     * Number of ghost node layers of each spatial dimension.
     * Note that the global boundary accounts for one ghost layer,
     * except in the periodic boundary condition, in which the
     * global boundary should be treated as normal node layers.
     */
    for (int s = 0; s < DIMS; ++s) {
        part->ng[s] = part->gl - 1;
    }
    /* adjust according to dimension collapse */
    if (OPTSPLIT == model->multidim) {
        switch (part->collapse) {
            case COLLAPSEN:
                break;
            case COLLAPSEX:
                part->ng[X] = 0;
                break;
            case COLLAPSEY:
                part->ng[Y] = 0;
                break;
            case COLLAPSEZ:
                part->ng[Z] = 0;
                break;
            case COLLAPSEXY:
                part->ng[X] = 0;
                part->ng[Y] = 0;
                break;
            case COLLAPSEXZ:
                part->ng[X] = 0;
                part->ng[Z] = 0;
                break;
            case COLLAPSEYZ:
                part->ng[Y] = 0;
                part->ng[Z] = 0;
                break;
            default:
                break;
        }
    }
    /* adjust for periodic boundary conditions */
    for (int s = 0, p = PWB; s < DIMS; ++s, p = p + 2) {
        if (PERIODIC == part->typeBC[p]) {
            part->ng[s] = part->gl;
        }
    }
    /* mesh and node number on each spatial dimension */
    for (int s = 0; s < DIMS; ++s) {
        /* ensure at least two inner cells per dimension */
        part->m[s] = MaxInt(part->m[s], 2);
        /* total number of nodes (including ghost nodes) */
        part->n[s] = part->m[s] + 1 + 2 * part->ng[s];
    }
    return;
}
/*
 * This function computes and initializes the values of necessary
 * parameters, and performs the normalization of some parameters.
 */
static void InitializeParameters(Time *time, Space *space, Model *model)
{
    Partition *const part = &(space->part);
    Geometry *const geo = &(space->geo);
    /* space */
    for (int s = 0; s < DIMS; ++s) {
        part->domain[s][MAX] = part->domain[s][MAX] / model->refL;
        part->domain[s][MIN] = part->domain[s][MIN] / model->refL;
        part->d[s] = (part->domain[s][MAX] - part->domain[s][MIN]) / (Real)(part->m[s]);
        part->dd[s] = 1.0 / part->d[s];
    }
    part->tinyL = 1.0e-6 * MinReal(part->d[X], MinReal(part->d[Y], part->d[Z]));
    part->tinyL = part->tinyL * part->tinyL; /* distance square based comparison */
    /* time */
    time->end = time->end * model->refV / model->refL;
    if (0 >= time->stepN) {
        time->stepN = INT_MAX;
    }
    time->dataC = time->restart;
    for (int n = 0; n < NPROBE; ++n) {
        if (0 >= time->dataN[n]) {
            time->dataN[n] = 0;
        }
        if (0 >= time->dataW[n]) {
            time->dataW[n] = INT_MAX;
        }
    }
    /* geometry */
    if (0 >= geo->sphN) {
        geo->sphN = 0;
    }
    if (0 >= geo->stlN) {
        geo->stlN = 0;
    }
    geo->totN = geo->sphN + geo->stlN;
    /* model */
    if (0 >= model->ibmLayer) {
        model->ibmLayer = INT_MAX;
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
     * Now replace some parameters with general forms that are valid
     * for both dimensional and nondimensional N-S equations, since
     * dimensional forms can be seen as normalized by reference 1.
     */
    model->gasR = 1.0 / (model->gamma * model->refMa * model->refMa);
    model->cv = model->gasR / (model->gamma - 1.0);
    return;
}
/* a good practice: end file with a newline */

