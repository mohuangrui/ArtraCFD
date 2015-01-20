/****************************************************************************
 * Preprocess                                                               *
 * Last-modified: 17 Jan 2015 12:49:44 PM
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
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ProgramMemoryRelease(Field *, Flux *, Space *, Particle *, Partition *);
static int FinalInformation(void);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * this is the overall postprocessing function
 */
int Postprocess(Field *field, Flux *flux, Space *space, Particle *particle,
        Partition *part)
{
    ProgramMemoryRelease(field, flux, space, particle, part);
    FinalInformation();
    return 0;
}
/*
 * this function together with some subfuctions realize the dynamic
 * memory release for each global pointer
 */
static int ProgramMemoryRelease(Field *field, Flux *flux, Space *space,
        Particle *particle, Partition *part)
{
    ShowInformation("Releasing memory back to system...");
    /* field variable related */
    RetrieveStorage(field->U);
    RetrieveStorage(field->Un);
    RetrieveStorage(field->Um);
    RetrieveStorage(field->Uo);
    /* flux related */
    RetrieveStorage(flux->Fx);
    RetrieveStorage(flux->Fy);
    RetrieveStorage(flux->Fz);
    RetrieveStorage(flux->Gx);
    RetrieveStorage(flux->Gy);
    RetrieveStorage(flux->Gz);
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
 * this function shows some final information
 */
static int FinalInformation(void)
{
    ShowInformation("Computing finished, successfully exit!");
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

