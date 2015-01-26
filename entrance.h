/****************************************************************************
 * Header File                                                              *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 ****************************************************************************/
/****************************************************************************
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_ENTRANCE_H_ /* if this is the first definition */
#define ARTRACFD_ENTRANCE_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Program entrance
 *
 * Parameter
 *      argc -- main arguments
 *      argv -- main arguments
 * Function
 *      Call a series of functions to handle the entering of the program.
 * Returns
 *      0 -- successful
 */
int ProgramEntrance(int argc, char *argv[], Command *);
#endif
/* a good practice: end file with a newline */

