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
#include "casedataloader.h"
#include "cfdparameters.h"
#include "domainpartition.h"
#include "geometryloader.h"
#include "gcibm.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ProgramMemoryAllocate(Field *, Space *);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * This is the overall preprocessing function
 */
int Preprocess(Field *field, Space *space, Particle *particle, Time *time, 
        Partition *part, Flow *flow)
{
    LoadCaseSettingData(space, time, flow, part);
    ComputeCFDParameters(space, time, flow);
    DomainPartition(part, space);
    ProgramMemoryAllocate(field, space);
    LoadGeometryData(particle, time);
    InitializeDomainGeometry(space);
    ComputeDomainGeometryGCIBM(space, particle, part);
    return 0;
}
/*
 * This function together with some subfuctions realize the dynamic
 * memory allocation for each global pointer. The storage retrieving
 * need to be done in the postprocessor.
 */
static int ProgramMemoryAllocate(Field *field, Space *space)
{
    ShowInformation("Allocating memory for program...");
    /*
     * Conservative flow variables: rho, rho_u, rho_v, rho_w, rho_eT
     */
    int idxMax = space->nMax * 5;
    field->U = AssignStorage(idxMax, "Real");
    field->Un = AssignStorage(idxMax, "Real");
    field->Um = AssignStorage(idxMax, "Real");
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

