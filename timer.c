/****************************************************************************
 * Program Performance Timer                                                *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines functions for timing program running time             *
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
double TockTime(Timer *timer) {
    Timer now; /* store current time at now */
    gettimeofday(&now, NULL);
    return (double)(now.tv_sec - timer->tv_sec) + 
        ((double)(now.tv_usec - timer->tv_usec) / 1000000.0);
}
/* a good practice: end file with a newline */

