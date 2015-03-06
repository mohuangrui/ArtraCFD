/****************************************************************************
 * Preprocess                                                               *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file functions as a postprocessor                                 *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "postprocess.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ProgramMemoryRelease(Field *, Space *, Particle *);
static int FinalInformation(void);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * This is the overall postprocessing function
 */
int Postprocess(Field *field, Space *space, Particle *particle)
{
    ProgramMemoryRelease(field, space, particle);
    FinalInformation();
    return 0;
}
/*
 * This function together with some subfuctions realize the dynamic
 * memory release for each global pointer
 */
static int ProgramMemoryRelease(Field *field, Space *space, Particle *particle)
{
    ShowInformation("Releasing memory back to system...");
    /* field variable related */
    RetrieveStorage(field->U);
    RetrieveStorage(field->Un);
    RetrieveStorage(field->Um);
    /* space related */
    RetrieveStorage(space->nodeFlag);
    /* particle related */
    RetrieveStorage(particle->headAddress);
    ShowInformation("Session End");
    return 0;
}
/*
 * This function shows some final information
 */
static int FinalInformation(void)
{
    ShowInformation("Computing finished, successfully exit!");
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

