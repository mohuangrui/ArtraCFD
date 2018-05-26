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
#include "timer.h"
#include <stdio.h> /* standard library for input and output */
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function record the time at the calling moment into the argument
 */
void TickTime(Timer *timer) {
    gettimeofday(timer, NULL);
}
/*
 * This function returns the time in seconds from now to time described by
 * the input argument
 */
double TockTime(const Timer *timer) {
    Timer now; /* store current time at now */
    gettimeofday(&now, NULL);
    return (double)(now.tv_sec - timer->tv_sec) + 
        ((double)(now.tv_usec - timer->tv_usec) / 1000000.0);
}
/* a good practice: end file with a newline */

