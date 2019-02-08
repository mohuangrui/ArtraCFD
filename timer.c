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
#include <stddef.h> /* standard library for macros */
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Record the time at the calling moment into the argument
 */
void TickTime(Timer *tm)
{
    gettimeofday(tm, NULL);
}
/*
 * Return the time in seconds from now to the time described by the argument
 */
double TockTime(const Timer *tm)
{
    Timer tc; /* store current time at now */
    gettimeofday(&tc, NULL);
    return (double)(tc.tv_sec - tm->tv_sec) +
        (double)(tc.tv_usec - tm->tv_usec) / 1000000.0;
}
/* a good practice: end file with a newline */

