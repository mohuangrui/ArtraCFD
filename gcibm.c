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
static int InitializeDomainGeometry(Space *, const Partition *);
static int LocateSolidGeometry(Space *, const Particle *, const Partition *);
static int IdentifyGhostNodes(Space *, const Partition *);
static int IdentifySolidNodeWithGhostNeighbours(Space *, const Particle *, const Partition *);
static int Min(const int x, const int y);
static int Max(const int x, const int y);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * These functions identify the type of each node: 
 * 1:                   boundary and exterior ghost node,
 * <= -offset:          interior solid node,
 * >=offset:            interior ghost node, 
 * 0:                   interior fluid node,
 * <= -offset - totalN: interior solid node with ghost neighbour.
 *
 * Procedures are:
 * -- initialize node flag of boundary and exterior nodes to boundary type,
 *    and inner nodes to fluid type;
 * -- identify all inner nodes that are in solid geometry as solid node;
 * -- identify ghost nodes according to the node type of its neighbours;
 * -- identify whether a solid node has ghost neighbours.
 *
 * It's necessary to difference boundary nodes and inner nodes because this
 * will make the identification of ghost nodes much easier in both 2D and
 * 3D situations.
 *
 * Moreover, whenever identifying a solid node or ghost node, store its
 * corresponding geometry information by linking the node flag to the
 * geometry ID, which will be accessed to calculate other informations. 
 * The rational is that don't store every information for each ghost node, but
 * only store necessary information. When need it, access and calculate it.
 */
int ComputeDomainGeometryGCIBM(Space *space, Particle *particle, const Partition *part)
{
    InitializeDomainGeometry(space, part);
    LocateSolidGeometry(space, particle, part);
    IdentifyGhostNodes(space, part);
    IdentifySolidNodeWithGhostNeighbours(space, particle, part);
    return 0;
}
static int InitializeDomainGeometry(Space *space, const Partition *part)
{
    /*
     * Set the vaule of offset to specify the range assignment for node type
     * identifier.
     */
    space->nodeFlagOffset = 10;
    /*
     * Initialize the entire domain to boundary type. Operation can be achieved
     * by a single loop since all data are stored by linear arrays.
     */
    int idx = 0; /* linear array index math variable */
    for (idx = 0; idx < space->nMax; ++idx) {
        space->nodeFlag[idx] = 1;
    }
    /*
     * Initialize inner nodes to fluid type.
     */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                space->nodeFlag[idx] = 0;
            }
        }
    }
    return 0;
}
/*
 * When locate solid nodes, there are two approaches available. One is search
 * over each node and verify each node regarding to all the particles; another
 * is search each particle and find all the nodes inside current particle.
 * The second method is adopted here for performance reason, although it's much
 * more complicated than the first one.
 */
static int LocateSolidGeometry(Space *space, const Particle *particle, const Partition *part)
{
    int idx = 0; /* linear array index math variable */
    /* geometry computation */
    Real distance = 0.0;
    Real distX = 0.0;
    Real distY = 0.0;
    Real distZ = 0.0;
    Real radius = 0.0;
    const Real *ptk = particle->headAddress;
    const int offset = space->nodeFlagOffset;
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = ptk + geoCount * particle->entryN; /* point to storage of current particle */
        const int iCenter = (int)((ptk[0] - space->xMin) * space->ddx) + space->ng;
        const int jCenter = (int)((ptk[1] - space->yMin) * space->ddy) + space->ng;
        const int kCenter = (int)((ptk[2] - space->zMin) * space->ddz) + space->ng;
        const Real safetyCoe = 1.5; /* zoom the search range */
        const int iRange = (int)(safetyCoe * ptk[3] * space->ddx);
        const int jRange = (int)(safetyCoe * ptk[3] * space->ddy);
        const int kRange = (int)(safetyCoe * ptk[3] * space->ddz);
        const int kSub = Max(kCenter - kRange, part->kSub[0]);
        const int kSup = Min(kCenter + kRange + 1, part->kSup[0]);
        const int jSub = Max(jCenter - jRange, part->jSub[0]);
        const int jSup = Min(jCenter + jRange + 1, part->jSup[0]);
        const int iSub = Max(iCenter - iRange, part->iSub[0]);
        const int iSup = Min(iCenter + iRange + 1, part->iSup[0]);
        for (int k = kSub; k < kSup; ++k) {
            for (int j = jSub; j < jSup; ++j) {
                for (int i = iSub; i < iSup; ++i) {
                    idx = (k * space->jMax + j) * space->iMax + i;
                    distX = space->xMin + (i - space->ng) * space->dx - ptk[0];
                    distY = space->yMin + (j - space->ng) * space->dy - ptk[1];
                    distZ = space->zMin + (k - space->ng) * space->dz - ptk[2];
                    distance = distX * distX + distY * distY + distZ * distZ - ptk[3] * ptk[3];
                    if (0 > distance) { /* in the solid geometry */
                        space->nodeFlag[idx] = -offset - geoCount; /* geometry are linked */
                    }
                }
            }
        }
    }
    return 0;
}
static int IdentifyGhostNodes(Space *space, const Partition *part)
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
    const int offset = space->nodeFlagOffset;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (-offset < space->nodeFlag[idx]) { /* it's not solid node */
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
                    space->nodeFlag[idx] = -space->nodeFlag[idx]; /* geometry information conserved */
                }
            }
        }
    }
    return 0;
}
static int IdentifySolidNodeWithGhostNeighbours(Space *space, const Particle *particle, const Partition *part)
{
    int idx = 0; /* linear array index math variable */
    int idxW = 0; /* index at West */
    int idxE = 0; /* index at East */
    int idxS = 0; /* index at South */
    int idxN = 0; /* index at North */
    int idxF = 0; /* index at Front */
    int idxB = 0; /* index at Back */
    const int offset = space->nodeFlagOffset;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (-offset < space->nodeFlag[idx]) { /* it's not solid node */
                    continue;
                }
                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxE = (k * space->jMax + j) * space->iMax + i + 1;
                idxS = (k * space->jMax + j - 1) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                idxF = ((k - 1) * space->jMax + j) * space->iMax + i;
                idxB = ((k + 1) * space->jMax + j) * space->iMax + i;
                if ((-offset >= space->nodeFlag[idxW]) && (-offset >= space->nodeFlag[idxE]) &&  
                        (-offset >= space->nodeFlag[idxS]) && (-offset >= space->nodeFlag[idxN]) && 
                        (-offset >= space->nodeFlag[idxF]) && (-offset >= space->nodeFlag[idxB])) {
                    continue; /* this solid node has no ghost neighbour */
                }
                /* exist at least one neighbour is ghost */
                space->nodeFlag[idx] = space->nodeFlag[idx] - particle->totalN; /* geometry information conserved */
            }
        }
    }
    return 0;
}
/*
 * Boundary condition for interior ghost nodes
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
static int Min(const int x, const int y)
{
    if (x < y) {
        return x;
    }
    return y;
}
static int Max(const int x, const int y)
{
    if (x > y) {
        return x;
    }
    return y;
}
/* a good practice: end file with a newline */

