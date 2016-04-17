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
     * Declare and initialize variables. Anything in C can be initialised
     * with = 0; this initialises numeric elements to zero and pointers null
     */    
    Space theSpace = {0};
    Time theTime = {0};
    Model theModel = {0};
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
    Preprocess(&theSpace, &theTime, &theModel);
    /*
     * Solve
     */
    Solve(&theSpace, &theTime, &theModel);
    /*
     * Postprocessing
     */
    Postprocess(&theSpace);
    /*
     * Successfully return
     */
    exit(EXIT_SUCCESS); /* exiting program */ 
}
/* a good practice: end file with a newline */

