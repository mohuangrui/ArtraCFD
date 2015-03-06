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
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Coordinate system: right-handed Cartesian system.
 * X-Y plane is the screen plane, X is horizontal from west to east, 
 * Y is vertical from south to north. Z axis is perpendicular to the
 * screen and points from front to back. 
 *
 * The entire space domain will be decomposed according to the specified number
 * of processors. Moreover, each partition domain will be slitted into 13 parts:
 * Interior Cells
 * Domain [West],[East],[South],[North],[Front],[Back] Boundary
 * [West],[East],[South],[North],[Front],[Back] Exterior Ghost
 */
int DomainPartition(Partition *part, const Space *space)
{
    ShowInformation("Domain partitioning...");
    /*
     * Give names of each part
     */
    const char partName[13][15] = {
        "Interior",
        "Domain West",       "Domain East", 
        "Domain South",      "Domain North", 
        "Domain Front",      "Domain Back",
        "West Ghost",        "East Ghost", 
        "South Ghost",       "North Ghost", 
        "Front Ghost",       "Back Ghost"};
    for (int partCount = 0; partCount < part->subN; ++partCount) {
        strncpy(part->name[partCount], partName[partCount], sizeof part->name[partCount]);
    }
    /*
     * Outward facing surface unit normal vector values of each inner part, the
     * introducing of surface normal vector can provide great advantange: every
     * surface can be handled uniformly if manipulations and calculations are
     * incorporated with surface normal vector. For example, (normalX, normalY,
     * normalZ) is the outward facing unit normal vector at surface point (i,
     * j, k), then its neighbour node (ih, jh, kh) is more inner than current
     * node if (ih-i)*normalX + (jh-j)*normalY + (kh-k)*normalZ < 0.
     */
    /* z component */
    part->normalZ[0] = 0;
    part->normalZ[1]  = 0;                            part->normalZ[2]  = 0;
    part->normalZ[3]  = 0;                            part->normalZ[4]  = 0;
    part->normalZ[5]  = -1;                           part->normalZ[6]  = 1;
    /* y componnet */
    part->normalY[0] = 0;
    part->normalY[1]  = 0;                            part->normalY[2]  = 0;
    part->normalY[3]  = -1;                           part->normalY[4]  = 1;
    part->normalY[5]  = 0;                            part->normalY[6]  = 0;
    /* x component */
    part->normalX[0] = 0;
    part->normalX[1]  = -1;                           part->normalX[2]  = 1;
    part->normalX[3]  = 0;                            part->normalX[4]  = 0;
    part->normalX[5]  = 0;                            part->normalX[6]  = 0;
    /*
     * Assign values to each index control array of each inner partition.
     *
     * Each index pair Sub[n] and Sup[n] defines the index ranges of the nth
     * inner partition. Remember that Sub is reachable value while Sup is 
     * unreachable value in for loop.
     *
     * Note that for each direction, its boundary nodes and exterior ghost 
     * nodes only need to extent out from the interior cells at that direction
     * and do not need to extent on other directions, that is, they form cross
     * like shapes in space without corner parts.
     */
    /* kSub */
    part->kSub[0]  = space->ng + 1;
    part->kSub[1]  = space->ng + 1;                   part->kSub[2]  = space->ng + 1;
    part->kSub[3]  = space->ng + 1;                   part->kSub[4]  = space->ng + 1;
    part->kSub[5]  = space->ng;                       part->kSub[6]  = space->nz + space->ng - 1;
    part->kSub[7]  = space->ng + 1;                   part->kSub[8]  = space->ng + 1;
    part->kSub[9]  = space->ng + 1;                   part->kSub[10] = space->ng + 1;
    part->kSub[11] = 0;                               part->kSub[12] = space->nz + space->ng;
    /* kSup */
    part->kSup[0]  = space->nz + space->ng - 1;
    part->kSup[1]  = space->nz + space->ng - 1;       part->kSup[2]  = space->nz + space->ng - 1;
    part->kSup[3]  = space->nz + space->ng - 1;       part->kSup[4]  = space->nz + space->ng - 1;
    part->kSup[5]  = space->ng + 1;                   part->kSup[6]  = space->nz + space->ng;
    part->kSup[7]  = space->nz + space->ng - 1;       part->kSup[8]  = space->nz + space->ng - 1;
    part->kSup[9]  = space->nz + space->ng - 1;       part->kSup[10] = space->nz + space->ng - 1;
    part->kSup[11] = space->ng;                       part->kSup[12] = space->nz + 2 * space->ng; 
    /* jSub */
    part->jSub[0]  = space->ng + 1;
    part->jSub[1]  = space->ng + 1;                   part->jSub[2]  = space->ng + 1;
    part->jSub[3]  = space->ng;                       part->jSub[4]  = space->ny + space->ng - 1;
    part->jSub[5]  = space->ng + 1;                   part->jSub[6]  = space->ng + 1;
    part->jSub[7]  = space->ng + 1;                   part->jSub[8]  = space->ng + 1;
    part->jSub[9]  = 0;                               part->jSub[10] = space->ny + space->ng;
    part->jSub[11] = space->ng + 1;                   part->jSub[12] = space->ng + 1;
    /* jSup */
    part->jSup[0]  = space->ny + space->ng - 1;
    part->jSup[1]  = space->ny + space->ng - 1;       part->jSup[2]  = space->ny + space->ng - 1;
    part->jSup[3]  = space->ng + 1;                   part->jSup[4]  = space->ny + space->ng;
    part->jSup[5]  = space->ny + space->ng - 1;       part->jSup[6]  = space->ny + space->ng - 1;
    part->jSup[7]  = space->ny + space->ng - 1;       part->jSup[8]  = space->ny + space->ng - 1;
    part->jSup[9]  = space->ng;                       part->jSup[10] = space->ny + 2 * space->ng; 
    part->jSup[11] = space->ny + space->ng - 1;       part->jSup[12] = space->ny + space->ng - 1;
    /* iSub */
    part->iSub[0]  = space->ng + 1;
    part->iSub[1]  = space->ng;                       part->iSub[2]  = space->nx + space->ng - 1;
    part->iSub[3]  = space->ng + 1;                   part->iSub[4]  = space->ng + 1;
    part->iSub[5]  = space->ng + 1;                   part->iSub[6]  = space->ng + 1;
    part->iSub[7]  = 0;                               part->iSub[8]  = space->nx + space->ng; 
    part->iSub[9]  = space->ng + 1;                   part->iSub[10] = space->ng + 1;
    part->iSub[11] = space->ng + 1;                   part->iSub[12] = space->ng + 1;
    /* iSup */
    part->iSup[0]  = space->nx + space->ng - 1;
    part->iSup[1]  = space->ng + 1;                   part->iSup[2]  = space->nx + space->ng;
    part->iSup[3]  = space->nx + space->ng - 1;       part->iSup[4]  = space->nx + space->ng - 1;
    part->iSup[5]  = space->nx + space->ng - 1;       part->iSup[6]  = space->nx + space->ng - 1;
    part->iSup[7]  = space->ng;                       part->iSup[8]  = space->nx + 2 * space->ng;  
    part->iSup[9]  = space->nx + space->ng - 1;       part->iSup[10] = space->nx + space->ng - 1;
    part->iSup[11] = space->nx + space->ng - 1;       part->iSup[12] = space->nx + space->ng - 1;
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

