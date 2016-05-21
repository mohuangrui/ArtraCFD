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
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_TIMER_H_ /* if this is the first definition */
#define ARTRACFD_TIMER_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include <sys/time.h> /* system time */
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef struct timeval Timer;
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Tick current time
 *
 * Function
 *      Records the time moment at calling into the input argument. This 
 *      function can be called before blocks of code to record current time
 *      which will be measured in microseconds. In short, tick the time and
 *      store in the timer.
 */
extern void TickTime(Timer *timer);
/*
 * Tock the timer
 *
 * Function
 *      Return the time interval between the time recorded by the timer and the
 *      moment at "now" (the calling time of this function) in the form of
 *      seconds. In short, tock the timer.
 */
extern double TockTime(const Timer *timer);
#endif
/* a good practice: end file with a newline */

