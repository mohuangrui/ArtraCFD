/****************************************************************************
 * Preprocess                                                               *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a functions as a preprocessor                        *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "preprocess.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <string.h> /* manipulating strings */
#include "casedataloader.h"
#include "cfdparameters.h"
#include "domainpartition.h"
#include "geometryloader.h"
#include "gcibm.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ProgramMemoryAllocate(Field *, Flux *, Space *);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * This is the overall preprocessing function
 */
int Preprocess(Field *field, Flux *flux, Space *space,
        Particle *particle, Time *time, Partition *part, Fluid *fluid, 
        Flow *flow, Reference *reference)
{
    LoadCaseSettingData(space, time, fluid, reference);
    ComputeCFDParameters(space, time, fluid, flow, reference);
    DomainPartition(part, space);
    ProgramMemoryAllocate(field, flux, space);
    LoadGeometryData(particle, time);
    InitializeDomainGeometryGCIBM(space, particle, part);
    return 0;
}
/*
 * This function together with some subfuctions realize the dynamic
 * memory allocation for each global pointer. The storage retrieving
 * need to be done in the postprocessor.
 */
static int ProgramMemoryAllocate(Field *field, Flux *flux, Space *space)
{
    ShowInformation("Allocating memory for program...");
    /*
     * Conservative flow variables: rho, rho_u, rho_v, rho_w, E
     */
    int dimU = 5; /* dimension of field variable U */
    int idxMax = dimU * space->nMax;
    field->U = AssignStorage(idxMax, "Real");
    field->Un = AssignStorage(idxMax, "Real");
    field->Um = AssignStorage(idxMax, "Real");
    flux->Fx = AssignStorage(idxMax, "Real");
    flux->Fy = AssignStorage(idxMax, "Real");
    flux->Fz = AssignStorage(idxMax, "Real");
    flux->Gx = AssignStorage(idxMax, "Real");
    flux->Gy = AssignStorage(idxMax, "Real");
    flux->Gz = AssignStorage(idxMax, "Real");
    /*
     * Primitive flow variables: rho, u, v, w, p, T
     */
    idxMax = (dimU + 1) * space->nMax;
    field->Uo = AssignStorage(idxMax, "Real");
    /*
     * Node type identifier
     */
    idxMax = space->nMax; /* max index */
    space->ghostFlag = AssignStorage(idxMax, "int");
    space->geoID = AssignStorage(idxMax, "int");
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

