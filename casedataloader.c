/****************************************************************************
 * Case Data Loader                                                         *
 * Last-modified: 19 Jan 2015 10:06:26 AM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a loader for the case settings                       *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "casedataloader.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int WriteVerifyData(const Space *, const Time *, const Fluid *,
        const Reference *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function load the case settings from the case file.
 * The key is to read and process file line by line. Use "*** begin" 
 * in the case file to identify and control the loading. 
 * The function scanf is notorious for its poor end-of-line handling,
 * Instead, use fgets to read a line of input and sscanf to process it.
 * Note: use a large enough number when using fgets to ensure reading
 * a whole line at a time. fgets will get the entire line including
 * the newline character (\n).
 * NOTE: if memory locations of input objects overlap, the behavior of
 * sscanf is undefined!
 * NOTE: sscanf can correctly handle any space in the target string as
 * well as in the format specifier, therefore, no need to process those
 * lines that will be processed by sscanf.
 */
int LoadCaseSettingData(Space *space, Time *time, Fluid *fluid, 
        Reference *reference)
{
    ShowInformation("Loading case setting data ...");
    FILE *filePointer = fopen("artracfd.case", "r");
    if (filePointer == NULL) {
        FatalError("Failed to open case data file");
    }
    /*
     * Read file line by line to get case setting data
     */
    char currentLine[200] = {'\0'}; /* store the current read line */
    int entryCount = 0; /* entry count */
    while (fgets(currentLine, sizeof currentLine, filePointer) != NULL) {
        CommandLineProcessor(currentLine); /* process current line */
        if (strncmp(currentLine, "space begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg, %lg, %lg", 
                    &(space->dx), &(space->dy), &(space->dz)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d, %d, %d", 
                    &(space->nx), &(space->ny), &(space->nz)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(space->ng)); 
            continue;
        }
        if (strncmp(currentLine, "time begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->restart)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg", &(time->totalTime)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg", &(time->numCFL)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->totalOutputTimes)); 
            continue;
        }
        if (strncmp(currentLine, "fluid begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg", &(fluid->density)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg", &(fluid->nu)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg", &(fluid->alpha)); 
            continue;
        }
        if (strncmp(currentLine, "reference begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg", &(reference->length)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg", &(reference->density)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg", &(reference->velocity)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%lg", &(reference->temperature)); 
            continue;
        }
    }
    fclose(filePointer); /* close current opened file */
    /*
     * Check missing information section in configuration
     */
    if (entryCount != 4) {
        FatalError("Missing or repeated necessary information section");
    }
    /*
     * Output loaded data for verification
     */
    WriteVerifyData(space, time, fluid, reference);
    ShowInformation("Session End");
    return 0;
}
/*
 * This function outputs the case setting data to a file for verification.
 */
static int WriteVerifyData(const Space *space, const Time *time, 
        const Fluid *fluid, const Reference *reference)
{
    ShowInformation("  Data outputted into artracfd.verify...");
    FILE *filePointer = fopen("artracfd.verify", "w");
    if (filePointer == NULL) {
        FatalError("Can not open file, write data failed");
    }
    /* output information to file */
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "# Input Conformation of ArtraCFD\n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "#       >> Space Domain << \n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "x, y, z length: %lg, %lg, %lg\n", 
            space->dx, space->dy, space->dz); 
    fprintf(filePointer, "x, y, z mesh number: %d, %d, %d\n", 
            space->nx, space->ny, space->nz); 
    fprintf(filePointer, "exterior ghost cell layers: %d\n", space->ng); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "#       >> Time Domain << \n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "restart flag: %d\n", time->restart); 
    fprintf(filePointer, "total evolution time: %lg\n", time->totalTime); 
    fprintf(filePointer, "CFL condition number: %lg\n", time->numCFL); 
    fprintf(filePointer, "exporting data times: %d\n", time->totalOutputTimes); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "#     >> Fluid Properties << \n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "fluid density: %lg\n", fluid->density); 
    fprintf(filePointer, "kinematic viscosity: %lg\n", fluid->nu); 
    fprintf(filePointer, "thermal diffusivity: %lg\n", fluid->alpha); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "#    >> Characteristic Values << \n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "length: %lg\n", reference->length); 
    fprintf(filePointer, "density: %lg\n", reference->density); 
    fprintf(filePointer, "velocity: %lg\n", reference->velocity); 
    fprintf(filePointer, "temperature: %lg\n", reference->temperature); 
    fprintf(filePointer, "#--------------------------------------\n\n"); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

