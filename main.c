/****************************************************************************
 * ArtraCFD Main Function                                                   *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide                                  *
 * - This is the main file of the CFD code, controls the overall program.   *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
#include "entrance.h"
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
    Field theField = { /* flow field variables */
        .U = NULL,
        .Un = NULL,
        .Um = NULL};
    Space theSpace = { /* space dimensions */
        .nx = 0,
        .ny = 0,
        .nz = 0,
        .ng = 0,
        .iMax = 0,
        .jMax = 0,
        .kMax = 0,
        .nMax = 0,
        .dx = 0.0,
        .dy = 0.0,
        .dz = 0.0,
        .xMin = 0.0,
        .yMin = 0.0,
        .zMin = 0.0,
        .xMax = 0.0,
        .yMax = 0.0,
        .zMax = 0.0,
        .ghostFlag = NULL,
        .geoID = NULL};
    Particle theParticle = { /* particle entities */
        .totalN = 0,
        .headAddress = NULL,
        .x = NULL,
        .y = NULL,
        .z = NULL,
        .r = NULL,
        .u = NULL,
        .v = NULL,
        .w = NULL};
    Time theTime = { /* time dimensions */
        .restart = 0,
        .totalTime = 0.0,
        .currentTime = 0.0,
        .dt = 0.0,
        .numCFL = 0.0,
        .totalStep = 0,
        .stepCount = 0,
        .totalOutputTimes = 0,
        .outputCount = 0};
    Flow theFlow = { /* flow parameters */
        .refMa = 0.0,
        .refMu = 0.0,
        .refPr = 0.0,
        .gamma = 0.0,
        .gasR = 0.0,
        .cv = 0.0,
        .refLength = 0.0,
        .refDensity = 0.0,
        .refVelocity = 0.0,
        .refTemperature = 0.0};
    Partition thePart = { /* domain partition control */
        .totalN = 1,
        .subN = 13,
        .kSub = {0},
        .kSup = {0},
        .jSub = {0},
        .jSup = {0},
        .iSub = {0},
        .iSup = {0},
        .normalZ = {0},
        .normalY = {0},
        .normalX = {0},
        .typeBC = {0},
        .valueBC = {{0.0}},
        .name = {" "},
        .typeIC = {0},
        .valueIC = {0.0}};
    Control theControl = { /* program overall control */
        .runMode = 'i',
        .processorN = 1};
    /*
     * Program Entrance
     */
    ProgramEntrance(argc, argv, &theControl);
    /*
     * Preprocessing
     */
    Preprocess(&theField, &theSpace, &theParticle, &theTime, &thePart, &theFlow);
    /*
     * Solve
     */
    Solve(&theField, &theSpace, &theParticle, &theTime, &thePart, &theFlow);
    /*
     * Postprocessing
     */
    Postprocess(&theField, &theSpace, &theParticle);
    /*
     * Successfully return
     */
    exit(EXIT_SUCCESS); /* exiting program */ 
}
/* a good practice: end file with a newline */

