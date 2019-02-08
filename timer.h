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
#ifndef ARTRACFD_TIMER_H_ /* if undefined */
#define ARTRACFD_TIMER_H_ /* set a unique marker */
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
 */
extern void TickTime(Timer *tm);
/*
 * Tock the timer
 */
extern double TockTime(const Timer *tm);
#endif
/* a good practice: end file with a newline */

