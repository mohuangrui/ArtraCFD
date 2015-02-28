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
static int ProgramMemoryRelease(Field *, Space *, Particle *, Partition *);
static int FinalInformation(void);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * This is the overall postprocessing function
 */
int Postprocess(Field *field, Space *space, Particle *particle,
        Partition *part)
{
    ProgramMemoryRelease(field, space, particle, part);
    FinalInformation();
    return 0;
}
/*
 * This function together with some subfuctions realize the dynamic
 * memory release for each global pointer
 */
static int ProgramMemoryRelease(Field *field, Space *space,
        Particle *particle, Partition *part)
{
    ShowInformation("Releasing memory back to system...");
    /* field variable related */
    RetrieveStorage(field->U);
    RetrieveStorage(field->Un);
    RetrieveStorage(field->Um);
    /* space related */
    RetrieveStorage(space->geoID);
    RetrieveStorage(space->ghostFlag);
    /* particle related */
    RetrieveStorage(particle->x);
    RetrieveStorage(particle->y);
    RetrieveStorage(particle->z);
    RetrieveStorage(particle->r);
    RetrieveStorage(particle->u);
    RetrieveStorage(particle->v);
    RetrieveStorage(particle->w);
    /* partition related */
    RetrieveStorage(part->nameHead);
    RetrieveStorage(part->idxHead);
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

