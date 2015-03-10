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
static int LocateSolidGeometry(Space *, const Particle *, const Partition *);
static int IdentifyGhostCells(Space *, const Partition *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * These functions identify the type of each node: 
 * boundary and exterior ghost cell (2), interior ghost cell (>= 10), 
 * interior solid cell (<= -10), interior fluid cells (0),
 * Procedures are: first, initialize nodeFlag of all nodes to boundary type;
 * then, identify all inner nodes as in solid geometry or fluid;
 * finally, identify ghost nodes according to its neighbours' node type;
 * moreover, whenever identifying a solid node or ghost node, store its
 * corresponding geometry information by linking the nodeFlag to the
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
     * Initialize the entire domain to type "2". Operation can be achieved
     * by a single loop since all data are stored by linear arrays.
     */
    for (int idx = 0; idx < space->nMax; ++idx) {
        space->nodeFlag[idx] = 2;
    }
    return 0;
}
int ComputeDomainGeometryGCIBM(Space *space, Particle *particle, const Partition *part)
{
    LocateSolidGeometry(space, particle, part);
    IdentifyGhostCells(space, part);
    return 0;
}
static int LocateSolidGeometry(Space *space, const Particle *particle, const Partition *part)
{
    int idx = 0; /* linear array index math variable */
    /* geometry computation */
    Real distance = 0.0;
    Real distX = 0.0;
    Real distY = 0.0;
    Real distZ = 0.0;
    Real radius = 0.0;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                space->nodeFlag[idx] = 0; /* reset to fluid */
                for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
                    radius = particle->r[geoCount];
                    distX = space->xMin + (i - space->ng) * space->dx - particle->x[geoCount];
                    distY = space->yMin + (j - space->ng) * space->dy - particle->y[geoCount];
                    distZ = space->zMin + (k - space->ng) * space->dz - particle->z[geoCount];
                    distance = distX * distX + distY * distY + distZ * distZ - radius * radius;
                    if (0 > distance) { /* in the solid geometry */
                        space->nodeFlag[idx] = -10 - geoCount; /* geoID are linked */
                    }
                }
            }
        }
    }
    return 0;
}
static int IdentifyGhostCells(Space *space, const Partition *part)
{
    int idx = 0; /* linear array index math variable */
    int idxW = 0; /* index at West */
    int idxE = 0; /* index at East */
    int idxS = 0; /* index at South */
    int idxN = 0; /* index at North */
    int idxF = 0; /* index at Front */
    int idxB = 0; /* index at Back */
    /* criteria */
    int flag = 0;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (-10 < space->nodeFlag[idx]) { /* it's not solid cell */
                    continue;
                }
                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxE = (k * space->jMax + j) * space->iMax + i + 1;
                idxS = (k * space->jMax + j - 1) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                idxF = ((k - 1) * space->jMax + j) * space->iMax + i;
                idxB = ((k + 1) * space->jMax + j) * space->iMax + i;
                flag = space->nodeFlag[idxW] * space->nodeFlag[idxE] * 
                    space->nodeFlag[idxS] * space->nodeFlag[idxN] * 
                    space->nodeFlag[idxF] * space->nodeFlag[idxB];
                if (0 == flag) { /* if exist one neighbour is fluid, then it's ghost */
                    space->nodeFlag[idx] = -space->nodeFlag[idx]; /* geoID information conserved */
                }
            }
        }
    }
    return 0;
}
/*
 * Boundary condition for interior ghost cells
 */
int BoundaryConditionGCIBM(Real *U, const Space *space, const Particle *particle, 
        const Partition *part)
{
    int idx = 0; /* linear array index math variable */
    int idxW = 0; /* index at West */
    int idxE = 0; /* index at East */
    int idxS = 0; /* index at South */
    int idxN = 0; /* index at North */
    int idxF = 0; /* index at Front */
    int idxB = 0; /* index at Back */
    int geoID = 0; /* geometry id */
    Real distToCenter = 0.0; /* distance from node to particle center */
    Real distToSurface = 0.0; /* distance from node to particle surface */
    Real distX = 0.0;
    Real distY = 0.0;
    Real distZ = 0.0;
    Real radius = 0.0;
    Real normalX = 0.0; /* x component of normal vector at surface */
    Real normalY = 0.0; /* y component of normal vector at surface */
    Real normalZ = 0.0; /* z component of normal vector at surface */
    int imageI = 0; /* node coordinates of the image point of the ghost */
    int imageJ = 0; /* node coordinates of the image point of the ghost */
    int imageK = 0; /* node coordinates of the image point of the ghost */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (10 > space->nodeFlag[idx]) { /* it's not a ghost */
                    continue;
                }
                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxE = (k * space->jMax + j) * space->iMax + i + 1;
                idxS = (k * space->jMax + j - 1) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                idxF = ((k - 1) * space->jMax + j) * space->iMax + i;
                idxB = ((k + 1) * space->jMax + j) * space->iMax + i;

                geoID = space->nodeFlag[idx] - 10; /* extract geoID from inner ghost node flag */
                radius = particle->r[geoID];
                distX = space->xMin + (i - space->ng) * space->dx - particle->x[geoID];
                distY = space->yMin + (j - space->ng) * space->dy - particle->y[geoID];
                distZ = space->zMin + (k - space->ng) * space->dz - particle->z[geoID];
                distToCenter = sqrt(distX * distX + distY * distY + distZ * distZ);
                normalX = distX / distToCenter;
                normalY = distY / distToCenter;
                normalZ = distZ / distToCenter;
                distToSurface = radius - distToCenter;
                imageI = i + (int)(2 * distToSurface * normalX * space->ddx);
                imageJ = j + (int)(2 * distToSurface * normalY * space->ddy);
                imageK = k + (int)(2 * distToSurface * normalZ * space->ddz);
            }
            /* integrate the forces on current particle */
        }
    }
    return 0;
}
/* a good practice: end file with a newline */

