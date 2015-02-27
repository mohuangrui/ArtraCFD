/****************************************************************************
 * Space Domain Meshing                                                     *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file use ghost cell immersed boundary method to handle complex    *
 *   geometries                                                             *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "gcibm.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int LocateSolidGeometry(Space *, Particle *, const Partition *);
static int IdentifyGhostCells(Space *, const Partition *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * These functions identify the type of each node: 
 * boundary and exterior ghost cell (2), interior ghost cell (1), 
 * interior solid cell (-1), interior fluid cells (0),
 * The procedures are: first, initialize ghostFlag of all nodes to type (2);
 * then, identify all inner nodes as in solid geometry (-1) or fluid (0);
 * finally, identify ghost nodes according to its neighbours' node type;
 * moreover, whenever identifying a solid node or ghost node, store its
 * corresponding geometry information by linking the node to the corresponding
 * geometry ID, which will be accessed to calculate other informations. 
 * The rational is that don't store every information for each ghost cell, but
 * only store necessary information. When need it, access and calculate it.
 */
/*
 * This function initializes the geometry for the entire domain.
 * It's separated because only need to be executed once as initialization.
 */
int InitializeDomainGeometry(Space *space)
{
    /*
     * Initialize the entire domain to type "2"
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* linear array index math variable */
    for (k = 0; k < space->kMax; ++k) {
        for (j = 0; j < space->jMax; ++j) {
            for (i = 0; i < space->iMax; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                space->ghostFlag[idx] = 2;
            }
        }
    }
    return 0;
}
int ComputeDomainGeometryGCIBM(Space *space, Particle *particle, const Partition *part)
{
    ShowInformation("  Computing domain geometry...");
    LocateSolidGeometry(space, particle, part);
    IdentifyGhostCells(space, part);
    return 0;
}
static int LocateSolidGeometry(Space *space, Particle *particle, const Partition *part)
{
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* linear array index math variable */
    /* geometry computation */
    int geoCount = 0; /* geometry objects count */
    Real distance = 0;
    Real distX = 0;
    Real distY = 0;
    Real distZ = 0;
    Real radius = 0;
    for (k = part->kSub[12]; k < part->kSup[12]; ++k) {
        for (j = part->jSub[12]; j < part->jSup[12]; ++j) {
            for (i = part->iSub[12]; i < part->iSup[12]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                space->ghostFlag[idx] = 0; /* reset to fluid */
                for (geoCount = 0; geoCount < particle->totalN; ++geoCount) {
                    radius = particle->r[geoCount];
                    distX = (i - space->ng) * space->dx - particle->x[geoCount];
                    distY = (j - space->ng) * space->dy - particle->y[geoCount];
                    distZ = (k - space->ng) * space->dz - particle->z[geoCount];
                    distance = sqrt(distX * distX + distY * distY + distZ * distZ) - radius;
                    if (distance < 0) { /* in the solid geometry */
                        space->ghostFlag[idx] = -1;
                        space->geoID[idx] = geoCount;
                    }
                }
            }
        }
    }
    return 0;
}
static int IdentifyGhostCells(Space *space, const Partition *part)
{
    /* indices */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* linear array index math variable */
    int idxW = 0; /* index at West */
    int idxE = 0; /* index at East */
    int idxS = 0; /* index at South */
    int idxN = 0; /* index at North */
    int idxF = 0; /* index at Front */
    int idxB = 0; /* index at Back */
    /* criteria */
    int flag = 0;
    for (k = part->kSub[12]; k < part->kSup[12]; ++k) {
        for (j = part->jSub[12]; j < part->jSup[12]; ++j) {
            for (i = part->iSub[12]; i < part->iSup[12]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (space->ghostFlag[idx] != -1) { /* it's not solid cell */
                    continue;
                }
                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxE = (k * space->jMax + j) * space->iMax + i + 1;
                idxS = (k * space->jMax + j - 1) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                idxF = ((k - 1) * space->jMax + j) * space->iMax + i;
                idxB = ((k + 1) * space->jMax + j) * space->iMax + i;
                flag = space->ghostFlag[idxW] * space->ghostFlag[idxE] * 
                    space->ghostFlag[idxS] * space->ghostFlag[idxN] * 
                    space->ghostFlag[idxF] * space->ghostFlag[idxB];
                if (flag == 0) { /* if exist one neighbour is fluid, then it's ghost */
                    space->ghostFlag[idx] = 1;
                }
            }
        }
    }
    return 0;
}
/*
 * Boundary condition for interior ghost cells
 */
int BoundaryConditionGCIBM(Field *field, const Space *space, const Particle *particle, 
        const Partition *part)
{
    /*
     * Decompose the field variable into each component.
     */
    Real *Un[5] = {
        field->Un + 0 * space->nMax,
        field->Un + 1 * space->nMax,
        field->Un + 2 * space->nMax,
        field->Un + 3 * space->nMax,
        field->Un + 4 * space->nMax};
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int dim = 0; /* dimension count of vectors */
    int idx = 0; /* linear array index math variable */
    int idxW = 0; /* index at West */
    int idxE = 0; /* index at East */
    int idxS = 0; /* index at South */
    int idxN = 0; /* index at North */
    int idxF = 0; /* index at Front */
    int idxB = 0; /* index at Back */
    Real distToCenter = 0; /* distance from node to particle center */
    Real distToSurface = 0; /* distance from node to particle surface */
    Real distX = 0;
    Real distY = 0;
    Real distZ = 0;
    Real radius = 0;
    Real normalX = 0; /* x component of normal vector at surface */
    Real normalY = 0; /* y component of normal vector at surface */
    Real normalZ = 0; /* z component of normal vector at surface */
    int imageI = 0; /* node coordinates of the image point of the ghost */
    int imageJ = 0; /* node coordinates of the image point of the ghost */
    int imageK = 0; /* node coordinates of the image point of the ghost */
    for (k = part->kSub[12]; k < part->kSup[12]; ++k) {
        for (j = part->jSub[12]; j < part->jSup[12]; ++j) {
            for (i = part->iSub[12]; i < part->iSup[12]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (space->ghostFlag[idx] != 1) { /* it's not a ghost */
                    continue;
                }
                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxE = (k * space->jMax + j) * space->iMax + i + 1;
                idxS = (k * space->jMax + j - 1) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                idxF = ((k - 1) * space->jMax + j) * space->iMax + i;
                idxB = ((k + 1) * space->jMax + j) * space->iMax + i;

                radius = particle->r[space->geoID[idx]];
                distX = (i - space->ng) * space->dx - particle->x[space->geoID[idx]];
                distY = (j - space->ng) * space->dy - particle->y[space->geoID[idx]];
                distZ = (k - space->ng) * space->dz - particle->z[space->geoID[idx]];
                distToCenter = sqrt(distX * distX + distY * distY + distZ * distZ);
                normalX = distX / distToCenter;
                normalY = distY / distToCenter;
                normalZ = distZ / distToCenter;
                distToSurface = radius - distToCenter;
                imageI = i + (int)(2 * distToSurface * normalX / space->dx);
                imageJ = j + (int)(2 * distToSurface * normalY / space->dy);
                imageK = k + (int)(2 * distToSurface * normalZ / space->dz);
            }
        }
    }
    return 0;
}
/* a good practice: end file with a newline */

