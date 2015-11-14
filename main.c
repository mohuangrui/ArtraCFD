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
#include "commons.h"
#include "program_entrance.h"
#include "preprocess.h"
#include "solve.h"
#include "postprocess.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <string.h> /* manipulating strings */
/****************************************************************************
 * The Main Function
 ****************************************************************************/
int main(int argc, char *argv[])
{
    /*
     * Declare and initialize variables
     */    
    Space theSpace = {
        .node = NULL,
        .m = {0},
        .n = {0},
        .totalN = 0,
        .ng = 0,
        .collapsed = 0,
        .d = {0.0},
        .dd = {0.0},
        .tinyL = 0.0,
        .domain = {{0.0}}
    };
    Time theTime = {
        .restart = 0,
        .stepN = 0,
        .stepCount = 0,
        .outputN = 0,
        .outputCount = 0,
        .dataStreamer = 0,
        .outputNProbe = 0,
        .tallyProbe = 0,
        .end = 0.0,
        .now = 0.0,
        .dt = 0.0,
        .numCFL = 0.0,
        .probe = {{0.0}}
    };
    Model theModel = {
        .scheme = 0,
        .averager = 0,
        .splitter = 0,
        .fsi = 0,
        .layers = 0,
        .fluid = 0,
        .refMa = 0.0,
        .refMu = 0.0,
        .gamma = 0.0,
        .gasR = 0.0,
        .cv = 0.0,
        .refLength = 0.0,
        .refDensity = 0.0,
        .refVelocity = 0.0,
        .refTemperature = 0.0
    };
    Partition thePart = {
        .n = {{{0}}},
        .normal = {{0}},
        .typeBC = {0},
        .tallyIC = 0,
        .typeIC = {0},
        .valueBC = {{0.0}},
        .valueIC = {{0.0}}
    };
    Geometry theGeo = {
        .totalM = 0,
        .sphereM = 0,
        .stlM = 0,
        .list = NULL
    };
    Control theControl = {
        .runMode = 'i',
        .procN = 1
    };
    /*
     * Program Entrance
     */
    ProgramEntrance(argc, argv, &theControl);
    /*
     * Preprocessing
     */
    Preprocess(&theSpace, &theTime, &theModel, &thePart, &theGeo);
    /*
     * Solve
     */
    Solve(&theSpace, &theTime, &theModel, &thePart, &theGeo);
    /*
     * Postprocessing
     */
    Postprocess(&theSpace, &theGeo);
    /*
     * Successfully return
     */
    exit(EXIT_SUCCESS); /* exiting program */ 
}
/* a good practice: end file with a newline */

