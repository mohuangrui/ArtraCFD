/****************************************************************************
 * Functions for Domain Decomposition                                       *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a function for partitioning computational domain.    *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "domainpartition.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Coordinate system: right-handed Cartesian system.
 * X-Y plane is the screen plane, X is horizontal from west to east, 
 * Y is vertical from south to north. Z axis is perpendicular to the
 * screen and points from front to back. The origin locates at the 
 * west-south-front corner of the computational domain.
 * This function decompose the space domain and provides index ranges. 
 * The current whole domain for a 3D space is slitted into 13 parts:
 * [West],[East],[South],[North],[Front],[Back] Exterior Ghost
 * Domain [West],[East],[South],[North],[Front],[Back] Boundary
 * Interior Cells
 */
int DomainPartition(Partition *part, const Space *space)
{
    ShowInformation("Domain partitioning...");
    /*
     * Specify the total number of domain partitions
     */
    part->totalN = 13; /* total number of partitions */
    /*
     * Give names of each part
     */
    const char partName[13][25] = {
        "West Exterior Ghost",        "East Exterior Ghost", 
        "South Exterior Ghost",       "North Exterior Ghost", 
        "Front Exterior Ghost",       "Back Exterior Ghost",
        "Domain West Boundary",       "Domain East Boundary", 
        "Domain South Boundary",      "Domain North Boundary", 
        "Domain Front Boundary",      "Domain Back Boundary",
        "Interior Cells"
    };
    /*
     * Assign storages and values for name pointers according to 
     * the total number of partitions and name length. These assigned
     * storage need to be received in the postprocessor.
     */
    /* the storage space for names */
    part->nameLength = 50; /* length of names */
    int idxMax = part->totalN * part->nameLength;
    part->nameHead = AssignStorage(idxMax, "char");
    /* distribute storage and assign names */
    int partCount = 0; /* count */
    int idx = 0; /* index math */
    for (partCount = 0; partCount < part->totalN; ++partCount) {
        idx = partCount * part->nameLength;
        strncpy(part->nameHead + idx, partName[partCount], part->nameLength);
    }
    /*
     * Assign storages and values for index pointers according to 
     * the total number of partitions
     */
    int idxDim = 6; /* number of index pointers */
    idxMax = idxDim * part->totalN;
    part->idxHead = AssignStorage(idxMax, "int");
    /* distribute storage to each index pointer */
    part->kSub = part->idxHead + 0 * part->totalN;
    part->kSup = part->idxHead + 1 * part->totalN;
    part->jSub = part->idxHead + 2 * part->totalN;
    part->jSup = part->idxHead + 3 * part->totalN;
    part->iSub = part->idxHead + 4 * part->totalN;
    part->iSup = part->idxHead + 5 * part->totalN;
    /*
     * Assign values to each index control array.
     *
     * Each index pair Sub[n] and Sup[n] defines the index ranges of the nth 
     * partition. Recall that Sub is reachable value while Sup is unreachable 
     * value in for loop.
     *
     * Note that for each direction, its boundary nodes and exterior ghost 
     * nodes only need to extent out from the interior cells at that direction
     * and do not need to extent on other directions, that is, they form cross
     * like shapes in space without corner parts.
     */
    /* kSub */
    part->kSub[0] = space->ng + 1;                    part->kSub[1] = space->ng + 1;
    part->kSub[2] = space->ng + 1;                    part->kSub[3] = space->ng + 1;
    part->kSub[4] = 0;                                part->kSub[5] = space->nz + space->ng;
    part->kSub[6] = space->ng + 1;                    part->kSub[7] = space->ng + 1;
    part->kSub[8] = space->ng + 1;                    part->kSub[9] = space->ng + 1;
    part->kSub[10] = space->ng;                       part->kSub[11] = space->nz + space->ng - 1;
    part->kSub[12] = space->ng + 1;
    /* kSup */
    part->kSup[0] = space->nz + space->ng - 1;        part->kSup[1] = space->nz + space->ng - 1;
    part->kSup[2] = space->nz + space->ng - 1;        part->kSup[3] = space->nz + space->ng - 1;
    part->kSup[4] = space->ng;                        part->kSup[5] = space->nz + 2 * space->ng; 
    part->kSup[6] = space->nz + space->ng - 1;        part->kSup[7] = space->nz + space->ng - 1;
    part->kSup[8] = space->nz + space->ng - 1;        part->kSup[9] = space->nz + space->ng - 1;
    part->kSup[10] = space->ng + 1;                   part->kSup[11] = space->nz + space->ng;
    part->kSup[12] = space->nz + space->ng - 1;
    /* jSub */
    part->jSub[0] = space->ng + 1;                    part->jSub[1] = space->ng + 1;
    part->jSub[2] = 0;                                part->jSub[3] = space->ny + space->ng;
    part->jSub[4] = space->ng + 1;                    part->jSub[5] = space->ng + 1;
    part->jSub[6] = space->ng + 1;                    part->jSub[7] = space->ng + 1;
    part->jSub[8] = space->ng;                        part->jSub[9] = space->ny + space->ng - 1;
    part->jSub[10] = space->ng + 1;                   part->jSub[11] = space->ng + 1;
    part->jSub[12] = space->ng + 1;
    /* jSup */
    part->jSup[0] = space->ny + space->ng - 1;        part->jSup[1] = space->ny + space->ng - 1;
    part->jSup[2] = space->ng;                        part->jSup[3] = space->ny + 2 * space->ng; 
    part->jSup[4] = space->ny + space->ng - 1;        part->jSup[5] = space->ny + space->ng - 1;
    part->jSup[6] = space->ny + space->ng - 1;        part->jSup[7] = space->ny + space->ng - 1;
    part->jSup[8] = space->ng + 1;                    part->jSup[9] = space->ny + space->ng;
    part->jSup[10] = space->ny + space->ng - 1;       part->jSup[11] = space->ny + space->ng - 1;
    part->jSup[12] = space->ny + space->ng - 1;
    /* iSub */
    part->iSub[0] = 0;                                part->iSub[1] = space->nx + space->ng; 
    part->iSub[2] = space->ng + 1;                    part->iSub[3] = space->ng + 1;
    part->iSub[4] = space->ng + 1;                    part->iSub[5] = space->ng + 1;
    part->iSub[6] = space->ng;                        part->iSub[7] = space->nx + space->ng - 1;
    part->iSub[8] = space->ng + 1;                    part->iSub[9] = space->ng + 1;
    part->iSub[10] = space->ng + 1;                   part->iSub[11] = space->ng + 1;
    part->iSub[12] = space->ng + 1;
    /* iSup */
    part->iSup[0] = space->ng;                        part->iSup[1] = space->nx + 2 * space->ng;  
    part->iSup[2] = space->nx + space->ng - 1;        part->iSup[3] = space->nx + space->ng - 1;
    part->iSup[4] = space->nx + space->ng - 1;        part->iSup[5] = space->nx + space->ng - 1;
    part->iSup[6] = space->ng + 1;                    part->iSup[7] = space->nx + space->ng;
    part->iSup[8] = space->nx + space->ng - 1;        part->iSup[9] = space->nx + space->ng - 1;
    part->iSup[10] = space->nx + space->ng - 1;       part->iSup[11] = space->nx + space->ng - 1;
    part->iSup[12] = space->nx + space->ng - 1;
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

