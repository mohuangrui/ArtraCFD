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
static int IdentifySolidNodesAtNumericalBoundary(Space *, const Particle *, 
        const Partition *);
static int SearchFluidNodes(const int k, const int j, const int i, 
        const int h, const Space *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * These functions identify the type of each node: 
 *
 * Procedures are:
 * -- initialize node flag of boundary and exterior nodes to boundary type,
 *    and inner nodes to fluid type;
 * -- identify all inner nodes that are in solid geometry as solid node;
 * -- identify ghost nodes according to the node type of its neighbours;
 * -- identify whether a solid node required for numerical boundary.
 *
 * It's necessary to difference boundary nodes and inner nodes because this
 * will not make incorrect mark of ghost nodes at nearby of domain boundaries.
 * It's necessary to difference ghost node and interior solid node which is
 * also required for numerical boundary because the latter is only required for
 * a numerical scheme has an accuracy order higher than 2nd. Besides, the ghost
 * node is special because the computation of boundary forces are only based on
 * the first nearest solid layer of boundary, that is, ghost node.
 *
 * The identification process should proceed step by step to correctly handle
 * all these relationships and avoid interference between each step,
 * interference may happen if identification processes are crunched together.
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
    IdentifySolidNodesAtNumericalBoundary(space, particle, part);
    return 0;
}
static int InitializeDomainGeometry(Space *space, const Partition *part)
{
    /*
     * Set the value of offset to specify the range assignment for node type
     * identifier.
     */
    space->nodeFlagOffset = 10;
    /*
     * Initialize the entire domain to boundary type. Operation can be achieved
     * by a single loop since all data are stored by linear arrays.
     */
    int idx = 0; /* linear array index math variable */
    for (idx = 0; idx < space->nMax; ++idx) {
        space->nodeFlag[idx] = -1;
    }
    /*
     * Initialize inner nodes to fluid type.
     */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
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
 * Be cautious with the validity of any calculated index. It's extremely
 * necessary to adjust the index into the valid flow region or check the
 * validity of the index to avoid index exceed array bound limits and 
 * mysterious bugs.
 */
static int LocateSolidGeometry(Space *space, const Particle *particle, const Partition *part)
{
    int idx = 0; /* linear array index math variable */
    Real info[10] = {0.0}; /* store calculated geometry information */
    const int offset = space->nodeFlagOffset;
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        const Real *ptk = particle->headAddress + geoCount * particle->entryN;
        const int iCenter = ComputeI(ptk[0], space);
        const int jCenter = ComputeJ(ptk[1], space);
        const int kCenter = ComputeK(ptk[2], space);
        const Real safetyCoe = 1.2; /* zoom the search range */
        const int iRange = (int)(safetyCoe * ptk[3] * space->ddx);
        const int jRange = (int)(safetyCoe * ptk[3] * space->ddy);
        const int kRange = (int)(safetyCoe * ptk[3] * space->ddz);
        /* determine search range according to valid flow region */
        const int kSub = FlowRegionK(kCenter - kRange, part);
        const int kSup = FlowRegionK(kCenter + kRange, part) + 1;
        const int jSub = FlowRegionJ(jCenter - jRange, part);
        const int jSup = FlowRegionJ(jCenter + jRange, part) + 1;
        const int iSub = FlowRegionI(iCenter - iRange, part);
        const int iSup = FlowRegionI(iCenter + iRange, part) + 1;
        for (int k = kSub; k < kSup; ++k) {
            for (int j = jSub; j < jSup; ++j) {
                for (int i = iSub; i < iSup; ++i) {
                    CalculateGeometryInformation(info, k, j, i, geoCount, space, particle);
                    if (0 < info[4]) { /* in the solid geometry */
                        idx = IndexMath(k, j, i, space);
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
    const int offset = space->nodeFlagOffset;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (-offset < space->nodeFlag[idx]) { /* it's not solid node */
                    continue;
                }
                /* if exist one neighbour is fluid, then it's ghost */
                if (0 == SearchFluidNodes(k, j, i, 1, space)) {
                    space->nodeFlag[idx] = -space->nodeFlag[idx]; /* geometry conserved */
                }
            }
        }
    }
    return 0;
}
static int IdentifySolidNodesAtNumericalBoundary(Space *space, 
        const Particle *particle, const Partition *part)
{
    int idx = 0; /* linear array index math variable */
    const int offset = space->nodeFlagOffset;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (-offset < space->nodeFlag[idx]) { /* it's not solid node */
                    continue;
                }
                for (int order = 2; order <= 2; ++order) {
                    if (0 == SearchFluidNodes(k, j, i, order, space)) {
                        space->nodeFlag[idx] = space->nodeFlag[idx] - particle->totalN;
                    }
                }
            }
        }
    }
    return 0;
}
/*
 * Search around current node and check whether current node has at
 * least one fluid node regarding to the specified coordinate range.
 */
static int SearchFluidNodes(const int k, const int j, const int i, 
        const int h, const Space *space)
{
    const int idxW = IndexMath(k, j, i - h, space);
    const int idxE = IndexMath(k, j, i + h, space);
    const int idxS = IndexMath(k, j - h, i, space);
    const int idxN = IndexMath(k, j + h, i, space);
    const int idxF = IndexMath(k - h, j, i, space);
    const int idxB = IndexMath(k + h, j, i, space);
    /* caution: enough ghost layers are quired to avoid illegal index */
    return (space->nodeFlag[idxW] * space->nodeFlag[idxE] * 
            space->nodeFlag[idxS] * space->nodeFlag[idxN] * 
            space->nodeFlag[idxF] * space->nodeFlag[idxB]);
}
/*
 * Boundary condition for interior ghost and solid nodes at numerical boundary.
 */
int BoundaryConditionGCIBM(Real *U, const Space *space, const Particle *particle, 
        const Partition *part, const Flow *flow)
{
    int idx = 0; /* linear array index math variable */
    Real Uo[6] = {0.0}; /* save reconstructed primitives */
    Real info[10] = {0.0}; /* store calculated geometry information */
    int geoID = 0; /* geometry id */
    const int offset = space->nodeFlagOffset;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space);
                if (-offset - particle->totalN >= space->nodeFlag[idx]) { /* solid */
                    geoID = space->nodeFlag[idx] + offset + particle->totalN;
                } else {
                    if (offset <= space->nodeFlag[idx]) { /* ghost */
                        geoID = space->nodeFlag[idx] - offset; /* extract geometry */
                    } else { /* not a numerical boundary node */
                        continue;
                    }
                }
                CalculateGeometryInformation(info, k, j, i, geoID, space, particle);
                /* obtain the spatial coordinates of the image point */
                const Real imageX = info[0] + 2 * info[4] * info[5];
                const Real imageY = info[1] + 2 * info[4] * info[6];
                const Real imageZ = info[2] + 2 * info[4] * info[7];
                Reconstruction(Uo, imageZ, imageY, imageX, 
                        ComputeK(imageZ, space), ComputeJ(imageY, space), ComputeI(imageX, space), 
                        U, space, flow);
                /*
                 * Apply no-slip wall boundary conditions to get primitive values at nodes
                 * in wall. That is, keep scalars and flip vectors after reflection.
                 * variable phi_ghost = 2 * phi_o - phi_image 
                 */
                Uo[1] = -Uo[1];
                Uo[2] = -Uo[2];
                Uo[3] = -Uo[3];
                ConservativeByPrimitive(U, idx * space->dimU, Uo, flow);
            }
        }
    }
    return 0;
}
/*
 * Reconstruction of the values of primitive vector Uo for a spatial
 * point (z, y, x) based on the neighbours around node (k, j, i). 
 * The inversed distance approach is adopted here.
 */
int Reconstruction(Real Uo[], const Real z, const Real y, const Real x,
        const int k, const int j, const int i, const Real *U,
        const Space *space, const Flow *flow)
{
    int idxh = 0; /* linear array index math variable */
    Real Uoh[6] = {0.0}; /* primitive at target node */
    Real Uow[6] = {0.0}; /* weighted primitive values */
    const Real tiny = 0.001 * (space->dx + space->dy + space->dz);
    Real weight = 0.0; /* weight factor */
    Real weightSum = 0.0; /* weight normalizer */
    /* 
     * Search around the specified node to find required fluid nodes as
     * interpolation stencil.
     */
    for (int kh = -2; kh < 3; ++kh) {
        for (int jh = -2; jh < 3; ++jh) {
            for (int ih = -2; ih < 3; ++ih) {
                idxh = IndexMath(k + kh, j + jh, i + ih, space);
                if ((0 > idxh) || (space->nMax <= idxh)) { /* check index */
                    continue; /* illegal index */
                }
                if (0 != space->nodeFlag[idxh]) { /* it's not a fluid node */
                    continue;
                }
                const Real distX = ComputeX(i + ih, space) - x;
                const Real distY = ComputeY(j + jh, space) - y;
                const Real distZ = ComputeZ(k + kh, space) - z;
                weight = distX * distX + distY * distY + distZ * distZ;
                if (tiny > weight) { /* avoid overflow of too small distance */
                    weight = tiny;
                }
                weight = 1 / weight;
                weightSum = weightSum + weight; /* accumulate normalizer */
                PrimitiveByConservative(Uoh, idxh * space->dimU, U, flow);
                for (int dim = 0; dim < space->dimU; ++dim) {
                    Uow[dim] = Uow[dim] + Uoh[dim] * weight;
                }
            }
        }
    }
    /*
     * Obtain the reconstructed primitive vector.
     */
    for (int dim = 0; dim < space->dimU; ++dim) {
        Uo[dim] = Uow[dim] / weightSum;
    }
    return 0;
}
int CalculateGeometryInformation(Real info[], const int k, const int j, const int i, 
        const int geoID, const Space *space, const Particle *particle)
{
    /* point to storage of current particle */
    const Real *ptk = particle->headAddress + geoID * particle->entryN;
    /* x, y, z coordinates */
    info[0] = ComputeX(i, space);
    info[1] = ComputeY(j, space);
    info[2] = ComputeZ(k, space);
    /* temporary store the x, y, z distance to particle center */
    info[5] = info[0] - ptk[0];
    info[6] = info[1] - ptk[1];
    info[7] = info[2] - ptk[2];
    /* distance to center */
    info[3] = sqrt(info[5] * info[5] + info[6] * info[6] + info[7] * info[7]);
    /* distance to surface, positive means in the geometry */
    info[4] = ptk[3] - info[3];
    /* x, y, z normal vector components */
    info[5] = info[5] / info[3];
    info[6] = info[6] / info[3];
    info[7] = info[7] / info[3];
    return 0;
}
/* a good practice: end file with a newline */

