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
#include "geometrystream.h"
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
    LoadGeometryData(particle, space, time, flow);
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
    ShowInformation("Allocating memory...");
    /*
     * Conservative flow variables.
     * Tips: the storage space of U is best between Un and Uswap.
     */
    int idxMax = space->nMax * DIMU;
    field->Un = AssignStorage(idxMax, "Real");
    field->U = AssignStorage(idxMax, "Real");
    field->Uswap = AssignStorage(idxMax, "Real");
    /*
     * Node type identifier
     */
    idxMax = space->nMax;
    space->nodeFlag = AssignStorage(idxMax, "int");
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

