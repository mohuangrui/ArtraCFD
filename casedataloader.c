/****************************************************************************
 * Case Data Loader                                                         *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a loader for the case setting file artracfd.case     *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "casedataloader.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ReadCaseSettingData(Space *, Time *, Flow *);
static int WriteVerifyData(const Space *, const Time *, const Flow *);
static int CheckCaseSettingData(const Space *, const Time *, const Flow *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function load the case settings from the case file.
 */
int LoadCaseSettingData(Space *space, Time *time, Flow *flow)
{
    ShowInformation("Loading case setting data ...");
    ReadCaseSettingData(space, time, flow);
    WriteVerifyData(space, time, flow);
    CheckCaseSettingData(space, time, flow);
    ShowInformation("Session End");
    return 0;
}
/*
 * This function read the case settings from the case file.
 * The key is to read and process file line by line. Use "*** begin" 
 * in the case file to identify and control the reading. 
 * The function scanf is notorious for its poor end-of-line handling.
 * Instead, use fgets to read a line of input and sscanf to process it.
 * Note: use a large enough number when using fgets to ensure reading
 * a whole line at a time. fgets will get the entire line including
 * the newline character (\n).
 * NOTE: if memory locations of input objects overlap, the behavior of
 * sscanf is undefined!
 * NOTE: sscanf can correctly handle any space in the target string as
 * well as in the format specifier, therefore, no need to process those
 * lines that will be processed by sscanf.
 * Footnote: In fprintf(), the rvalue type promotions are expected. %f and 
 * %g actually correspond to parameters of type double. Thus in fprintf()
 * there is no difference between %f and %lf, or between %g and %lg. However, 
 * in sscanf() what is passed is a pointer to the variable so no rvalue type 
 * promotions occur or are expected. Thus %f and %lf are quite different in
 * sscanf, but the same in fprintf. Consequently, we need to use %g for 
 * double in fprintf and %lg for double in sscanf. It doesn't matter which
 * you use for fprintf because the fprintf library function treats them as
 * synonymous, but it's crucial to get it right for sscanf. 
 */
static int ReadCaseSettingData(Space *space, Time *time, Flow *flow)
{
    FILE *filePointer = fopen("artracfd.case", "r");
    if (filePointer == NULL) {
        FatalError("failed to open case data file: artracfd.case...");
    }
    /*
     * Read file line by line to get case setting data
     */
    char currentLine[200] = {'\0'}; /* store the current read line */
    int entryCount = 0; /* entry count */
    /* set format specifier according to the type of Real */
    char formatI[5] = "%lg"; /* default is double type */
    char formatIII[15] = "%lg, %lg, %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
        strncpy(formatIII, "%g, %g, %g", sizeof formatIII); /* float type */
    }
    while (fgets(currentLine, sizeof currentLine, filePointer) != NULL) {
        CommandLineProcessor(currentLine); /* process current line */
        if (strncmp(currentLine, "space begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, 
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
            sscanf(currentLine, formatI, &(time->totalTime)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(time->numCFL)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->totalOutputTimes)); 
            continue;
        }
        if (strncmp(currentLine, "fluid begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refPr)); 
            continue;
        }
        if (strncmp(currentLine, "reference begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refLength)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refDensity)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refVelocity)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refTemperature)); 
            continue;
        }
    }
    fclose(filePointer); /* close current opened file */
    /*
     * Check missing information section in configuration
     */
    if (entryCount != 4) {
        FatalError("missing or repeated necessary information section");
    }
    return 0;
}
/*
 * This function outputs the case setting data to a file for verification.
 */
static int WriteVerifyData(const Space *space, const Time *time, 
        const Flow *flow)
{
    ShowInformation("  Data outputted into artracfd.verify...");
    FILE *filePointer = fopen("artracfd.verify", "w");
    if (filePointer == NULL) {
        FatalError("failed to write data to file: artracfd.verify");
    }
    /* output information to file */
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "# Input Conformation of ArtraCFD\n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "#       >> Space Domain << \n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "x, y, z length: %.6g, %.6g, %.6g\n", 
            space->dx, space->dy, space->dz); 
    fprintf(filePointer, "x, y, z mesh number: %d, %d, %d\n", 
            space->nx, space->ny, space->nz); 
    fprintf(filePointer, "exterior ghost cell layers: %d\n", space->ng); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "#       >> Time Domain << \n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "restart flag: %d\n", time->restart); 
    fprintf(filePointer, "total evolution time: %.6g\n", time->totalTime); 
    fprintf(filePointer, "CFL condition number: %.6g\n", time->numCFL); 
    fprintf(filePointer, "exporting data times: %d\n", time->totalOutputTimes); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "#  >> Fluid and Flow Properties << \n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "Prandtl number: %.6g\n", flow->refPr); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "#    >> Characteristic Values << \n"); 
    fprintf(filePointer, "#--------------------------------------\n"); 
    fprintf(filePointer, "length: %.6g\n", flow->refLength); 
    fprintf(filePointer, "density: %.6g\n", flow->refDensity); 
    fprintf(filePointer, "velocity: %.6g\n", flow->refVelocity); 
    fprintf(filePointer, "temperature: %.6g\n", flow->refTemperature); 
    fprintf(filePointer, "#--------------------------------------\n\n"); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * This function do some parameter checking
 */
static int CheckCaseSettingData(const Space *space, const Time *time, 
        const Flow *flow)
{
    ShowInformation("  Preliminary case data checking ...");
    /* space */
    if ((space->dz < 0) || (space->dy < 0) || (space->dx < 0)) {
        FatalError("negative length values in case settings");
    }
    if ((space->nz < 1) || (space->ny < 1) || (space->nx < 1)
            || (space->ng < 1)) {
        FatalError("too small mesh values in case settings");
    }
    /* time */
    if ((time->restart < 0) || (time->restart > 1)|| (time->totalTime <= 0)
            || (time->numCFL <= 0) || (time->totalOutputTimes < 1)) {
        FatalError("wrong values in time section of case settings");
    }
    /* fluid and flow */
    if ((flow->refPr <= 0)) {
        FatalError("wrong values in fluid and flow section of case settings");
    }
    /* reference */
    if ((flow->refLength <= 0) || (flow->refDensity <= 0) || 
            (flow->refVelocity <= 0) || flow->refTemperature <= 0) {
        FatalError("wrong values in reference section of case settings");
    }
    return 0;
}
/* a good practice: end file with a newline */

