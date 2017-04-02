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
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include "commons.h"
#include "program_entrance.h"
#include "preprocess.h"
#include "solve.h"
#include "postprocess.h"
/****************************************************************************
 * The Main Function
 ****************************************************************************/
int main(int argc, char *argv[])
{
    /*
     * Declare and initialize variables
     */    
    Control control = {
        .runMode = 'i',
        .procN = 1};
    Time time = {0};
    Space space = {0};
    Model model = {0};
    /*
     * Perform computation
     */
    ProgramEntrance(argc, argv, &control);
    Preprocess(&time, &space, &model);
    Solve(&time, &space, &model);
    Postprocess(&time, &space);
    exit(EXIT_SUCCESS); 
}
/* a good practice: end file with a newline */

