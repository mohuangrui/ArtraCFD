/****************************************************************************
 * Ensight Data Loader                                                      *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a function that load computed data                   *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "ensight.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int LoadEnsightCaseFile(EnsightSet *, Time *);
static int LoadEnsightVariableFile(Real *U, EnsightSet *,
        const Space *, const Partition *, const Flow *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Load necessary flow information from computed data
 */
int LoadComputedDataEnsight(Real *U, const Space *space, Time *time,
        const Partition *part, const Flow *flow)
{
    EnsightSet enSet = { /* initialize Ensight environment */
        .baseName = "restart", /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .stringData = {'\0'}, /* string data recorder */
    };
    LoadEnsightCaseFile(&enSet, time);
    LoadEnsightVariableFile(U, &enSet, space, part, flow);
    return 0;
}
static int LoadEnsightCaseFile(EnsightSet *enSet, Time *time)
{
    FILE *filePointer = NULL;
    /* current filename */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.case", enSet->baseName); 
    filePointer = fopen(enSet->fileName, "r");
    if (filePointer == NULL) {
        FatalError("failed to open restart case file: restart.case...");
    }
    /* read information from file */
    char currentLine[200] = {'\0'}; /* store current line */
    char garbage[100] = {'\0'}; /* store redundant information */
    /* set format specifier according to the type of Real */
    char formatV[25] = "%s %s %s %s %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatV, "%s %s %s %s %g", sizeof formatV); /* float type */
    }
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
    sscanf(currentLine, "%s %s %s %s %d", garbage, garbage,
            garbage, garbage, &(time->outputCount)); 
    /* get restart time */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, formatV, garbage, garbage,
            garbage, garbage, &(time->currentTime)); 
    /* get current step number */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%s %s %s %s %d", garbage, garbage,
            garbage, garbage, &(time->stepCount)); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int LoadEnsightVariableFile(Real *U, EnsightSet *enSet,
        const Space *space, const Partition *part, const Flow *flow)
{
    FILE *filePointer = NULL;
    int partCount = 0; /* part count starts from 0 */
    int partNum = 1; /* part number starts from 1 */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* linear array index math variable */
    EnsightReal data = 0; /* the ensight data format */
    int dim = 0; /* dimension count */
    const char nameSuffix[5][10] = {"den", "u", "v", "w", "pre"};
    /*
     * Define the primitive field variables.
     */
    Real rho = 0; 
    Real u = 0;
    Real v = 0;
    Real w = 0;
    for (dim = 0; dim < 5; ++dim) {
        snprintf(enSet->fileName, sizeof(EnsightString), "%s.%s", enSet->baseName,
                nameSuffix[dim]);
        filePointer = fopen(enSet->fileName, "rb");
        if (filePointer == NULL) {
            FatalError("failed to open restart data files: restart.***...");
        }
        fread(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        for (partCount = 0; partCount < part->totalN; ++partCount) {
            fread(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            fread(&partNum, sizeof(int), 1, filePointer);
            fread(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            for (k = part->kSub[partCount]; k < part->kSup[partCount]; ++k) {
                for (j = part->jSub[partCount]; j < part->jSup[partCount]; ++j) {
                    for (i = part->iSub[partCount]; i < part->iSup[partCount]; ++i) {
                        fread(&data, sizeof(EnsightReal), 1, filePointer);
                        idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                        switch (dim) {
                            case 0: /* rho */
                                U[idx] = data;
                                break;
                            case 1: /* u */
                                U[idx+1] = U[idx] * data;
                                break;
                            case 2: /* v */
                                U[idx+2] = U[idx] * data;
                                break;
                            case 3: /* w */
                                U[idx+3] = U[idx] * data;
                                break;
                            case 4: /* p */
                                rho = U[idx];
                                u = U[idx+1] / rho;
                                v = U[idx+2] / rho;
                                w = U[idx+3] / rho;
                                U[idx] = data / (flow->gamma - 1) + 0.5 * rho * (u * u + v * v + w * w);
                                break;
                            default:
                                break;
                        }
                    }
                }
            }
        }
        fclose(filePointer); /* close current opened file */
    }
    return 0;
}
/* a good practice: end file with a newline */

