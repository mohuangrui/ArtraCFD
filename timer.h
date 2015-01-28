/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
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
 * Parameter
 *      timer -- stores the time moment at calling
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
 * Parameter
 *      timer -- a timer which stores a past time moment.
 * Function
 *      Return the time interval between the time recorded by the timer and the
 *      moment at "now" (the calling time of this function) in the form of
 *      seconds. In short, tock the timer.
 */
extern double TockTime(Timer *timer);
#endif
/* a good practice: end file with a newline */

