/****************************************************************************
 * Flow Initialization                                                      *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a function that handles flow initialization.         *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "initialization.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "ensight.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int FirstRunInitializer(Field *, Space *, const Partition *);
static int RestartInitializer(Field *, Space *, Time *, const Partition *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function initializes the entire flow field. Initialization will be 
 * done differently determined by the restart status.
 */
int InitializeFlowField(Field *field, Space *space, const Particle *particle,
        Time *time, const Partition *part)
{
    ShowInformation("Initializing flow field...");
    if (time->restart == 0) { /* non restart */
        FirstRunInitializer(field, space, part);
        /* if this is a first run, output initial data */
        InitializeEnsightTransientCaseFile(time);
        WriteComputedDataEnsight(field->Uo, space, particle, time, part);
    } else {
        RestartInitializer(field, space, time,  part);
    }
    ShowInformation("Session End");
    return 0;
}
/*
 * The first run initialization will assign values to field variables.
 */
static int FirstRunInitializer(Field *field, Space *space, const Partition *part)
{
    ShowInformation("  Non-restart run initializing...");
    /*
     * Initialize flow field
     */
    return 0;
}
/*
 * If this is a restart run, then initialize flow field by reading field data
 * from restart files.
 */
static int RestartInitializer(Field *field, Space *space, Time *time,
        const Partition *part)
{
    ShowInformation("  Restart run initializing...");
    /*
     * Read the restart time and step information from restart.case file.
     */
    FILE *filePointer = NULL;
    const char restartBaseFileName[20] = "restart";
    char restartFileName[50] = {'\0'};
    snprintf(restartFileName, sizeof restartFileName, "%s.case", restartBaseFileName);
    filePointer = fopen(restartFileName, "r");
    if (filePointer == NULL) {
        FatalError("failed to open restart case file: restart.case...");
    }
    /* read information from file */
    char currentLine[200] = {'\0'}; /* store current line */
    char garbage[100] = {'\0'}; /* store redundant information */
    /* get rid of redundant lines */
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    fgets(currentLine, sizeof currentLine, filePointer);
    /* get restart order number */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%s %s %s %s %d\n", garbage, garbage,
            garbage, garbage, &(time->outputCount)); 
    /* get restart time */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%s %s %s %s %g\n", garbage, garbage,
            garbage, garbage, &(time->currentTime)); 
    /* get current step number */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%s %s %s %s %d\n", garbage, garbage,
            garbage, garbage, &(time->stepCount)); 
    fclose(filePointer); /* close current opened file */
    /*
     * Restore field data
     */
    int partCount = 0; /* part count starts from 0 */
    int partNum = 1; /* part number starts from 1 */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    /* linear array index math variables */
    int idx = 0; /* calculated index */
    const int stringLength = 80;
    const int dimU = 6; /* dimension of the primitive variable vector */
    int dimCount = 0; /* dimension count */
    const char nameSuffix[6][10] = {"den", "u", "v", "w", "pre", "tem"};
    for (dimCount = 0; dimCount < dimU; ++dimCount) {
        snprintf(restartFileName, sizeof restartFileName, "%s.%s",
                restartBaseFileName, nameSuffix[dimCount]);
        filePointer = fopen(restartFileName, "rb");
        if (filePointer == NULL) {
            FatalError("failed to open restart data files: restart.***...");
        }
        fread(currentLine, sizeof(char), stringLength, filePointer);
        for (partCount = 0; partCount < part->totalN; ++partCount) {
            fread(currentLine, sizeof(char), stringLength, filePointer);
            fread(&partNum, sizeof(int), 1, filePointer);
            fread(currentLine, sizeof(char), stringLength, filePointer);
            for (k = part->kSub[partCount]; k < part->kSup[partCount]; ++k) {
                for (j = part->jSub[partCount]; j < part->jSup[partCount]; ++j) {
                    for (i = part->iSub[partCount]; i < part->iSup[partCount]; ++i) {
                        idx = ((dimCount * space->kMax + k) * space->jMax + j) * space->iMax + i;
                        fread(field->Uo + idx, sizeof(float), 1, filePointer);
                    }
                }
            }
        }
        fclose(filePointer); /* close current opened file */
    }
    return 0;
}
/* a good practice: end file with a newline */

