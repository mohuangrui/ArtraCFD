/****************************************************************************
 * Parasight Data Loader                                                      *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a function that load computed data                   *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "parasight.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int LoadParasightCaseFile(ParasightSet *, Time *);
static int LoadParasightVariableFile(Real *U, ParasightSet *,
        const Space *, const Partition *, const Flow *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Load necessary flow information from computed data
 */
int LoadComputedDataParasight(Real *U, const Space *space, Time *time,
        const Partition *part, const Flow *flow)
{
    ParasightSet enSet = { /* initialize Parasight environment */
        .baseName = "restart", /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .stringData = {'\0'}, /* string data recorder */
    };
    LoadParasightCaseFile(&enSet, time);
    LoadParasightVariableFile(U, &enSet, space, part, flow);
    return 0;
}
static int LoadParasightCaseFile(ParasightSet *enSet, Time *time)
{
    FILE *filePointer = NULL;
    /* current filename */
    snprintf(enSet->fileName, sizeof(ParasightString), "%s.case", enSet->baseName); 
    filePointer = fopen(enSet->fileName, "r");
    if (NULL == filePointer) {
        FatalError("failed to open restart case file: restart.case...");
    }
    /* read information from file */
    char currentLine[200] = {'\0'}; /* store current line */
    /* set format specifier according to the type of Real */
    char format[25] = "%*s %*s %*s %*s %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(format, "%*s %*s %*s %*s %g", sizeof format); /* float type */
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
    sscanf(currentLine, "%*s %*s %*s %*s %d", &(time->outputCount)); 
    /* get restart time */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, format, &(time->currentTime)); 
    /* get current step number */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "%*s %*s %*s %*s %d", &(time->stepCount)); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
static int LoadParasightVariableFile(Real *U, ParasightSet *enSet,
        const Space *space, const Partition *part, const Flow *flow)
{
    FILE *filePointer = NULL;
    int idx = 0; /* linear array index math variable */
    ParasightReal data = 0.0; /* the Parasight data format */
    const char nameSuffix[5][10] = {"rho", "u", "v", "w", "p"};
    int partNum = 1;
    for (int dim = 0; dim < DIMU; ++dim) {
        snprintf(enSet->fileName, sizeof(ParasightString), "%s.%s", enSet->baseName, nameSuffix[dim]);
        filePointer = fopen(enSet->fileName, "rb");
        if (NULL == filePointer) {
            FatalError("failed to open restart data files...");
        }
        fread(enSet->stringData, sizeof(char), sizeof(ParasightString), filePointer);
        fread(enSet->stringData, sizeof(char), sizeof(ParasightString), filePointer);
        fread(&partNum, sizeof(int), 1, filePointer);
        fread(enSet->stringData, sizeof(char), sizeof(ParasightString), filePointer);
        for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
            for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
                for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                    fread(&data, sizeof(ParasightReal), 1, filePointer);
                    idx = IndexMath(k, j, i, space) * DIMU;
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
                            U[idx+4] = 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx] + data / flow->gammaMinusOne;
                            break;
                        default:
                            break;
                    }
                }
            }
        }
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

