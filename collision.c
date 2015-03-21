/****************************************************************************
 * Particle Collision Model                                                 *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "collision.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include "gcibm.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int SurfaceForceIntegration(const Real *U, const Space *space,
        Particle *particle, const Partition *part, const Flow *flow);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ParticleSpatialEvolution(Real *U, const Real dt, Space *space, Particle *particle, 
        const Partition *part, const Flow *flow)
{
    /*
     * Compute the forces acting on particles.
     */
    SurfaceForceIntegration(U, space, particle, part, flow);
    /*
     * Update particle velocity and position
     */
    int idx = 0; /* linear array index math variable */
    int geoID = 0; /* geometry id */
    Real *ptk = NULL;
    Real distance = 0.0;
    Real distX = 0.0;
    Real distY = 0.0;
    Real distZ = 0.0;
    Real mass = 0.0;
    const Real pi = acos(-1);
    const int offset = space->nodeFlagOffset;
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = particle->headAddress + geoCount * particle->entryN; /* point to storage of current particle */
        if (0 == space->dx * space->dy * space->dz) {
            mass = ptk[4] * pi * ptk[3] * ptk[3];
        } else {
            mass = ptk[4] * (4.0 / 3.0) * pi * ptk[3] * ptk[3] * ptk[3];
        }
        /* velocity */
        ptk[5] = ptk[5] + dt * ptk[8] / mass;
        ptk[6] = ptk[6] + dt * ptk[9] / mass;
        ptk[7] = ptk[7] + dt * ptk[10] / mass;
        /* spatial position */
        ptk[0] = ptk[0] + ptk[5] * dt;
        ptk[1] = ptk[1] + ptk[6] * dt;
        ptk[2] = ptk[2] + ptk[7] * dt;
    }
    /*
     * After the spatial positions of particles updated, some inner nodes fall
     * out of the solid region, these nodes need to be identified and the flow
     * values of these nodes need to be interpolated. Since the speed of
     * particles will not greater than the speed of fluid, the spatial motions
     * of particles will not exceed one grid per time step based on the CFL
     * condition of flow. Therefore, only ghost nodes need to be verified.
     */
    /* build an adequate search path */
    const int path[26][3] = { /* n paths for i, j, k */
        {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, 0}, {0, -1, 0}, {0, 0, -1},
        {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {-1, -1, 0}, {-1, 0, -1}, {0, -1, -1},
        {1, 1, 1}, {-1, -1, -1}, {1, -1, 0}, {1, 0, -1}, {0, 1, -1}, {0, -1, 1},
        {-1, 1, 0}, {-1, 0, 1}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}, {1, -1, -1},
        {-1, 1, -1}, {-1, -1, 1}};
    const int stencilN = 2; /* number of stencils for interpolation */
    int idxh = 0; /* index variable */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (offset > space->nodeFlag[idx]) { /* it's not a ghost */
                    continue;
                }
                geoID = space->nodeFlag[idx] - offset; /* extract geometry number from inner ghost node flag */
                ptk = particle->headAddress + geoID * particle->entryN; /* point to storage of current particle */
                distX = space->xMin + (i - space->ng) * space->dx - ptk[0];
                distY = space->yMin + (j - space->ng) * space->dy - ptk[1];
                distZ = space->zMin + (k - space->ng) * space->dz - ptk[2];
                distance = distX * distX + distY * distY + distZ * distZ - ptk[3] * ptk[3];
                if (0 > distance) { /* still in the solid geometry */
                    continue;
                }
                /* 
                 * Search around the image node to find required fluid nodes as
                 * interpolation stencil.
                 */
                Real Uo[5] = {0.0}; /* store weighted primitives */
                int tally = 0; /* number of current stencil */
                for (int loop = 0; (tally < stencilN) && (loop < 26); ++loop) {
                    const int ih = i + path[loop][0];
                    const int jh = j + path[loop][1];
                    const int kh = k + path[loop][2];
                    idxh = (kh * space->jMax + jh) * space->iMax + ih;
                    if (0 != space->nodeFlag[idxh]) { /* it's not a fluid node */
                        continue;
                    }
                    idxh = idxh * 5; /* switch to index field variable */
                    const Real rho_h = U[idxh+0];
                    const Real u_h = U[idxh+1] / rho_h;
                    const Real v_h = U[idxh+2] / rho_h;
                    const Real w_h = U[idxh+3] / rho_h;
                    const Real eT_h = U[idxh+4] / rho_h;
                    const Real p_h = (flow->gamma - 1.0) * rho_h * (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h));
                    Uo[0] = Uo[0] + (1 / stencilN) * rho_h;
                    Uo[1] = Uo[1] + (1 / stencilN) * u_h;
                    Uo[2] = Uo[2] + (1 / stencilN) * v_h;
                    Uo[3] = Uo[3] + (1 / stencilN) * w_h;
                    Uo[4] = Uo[4] + (1 / stencilN) * p_h;
                    ++tally; /* increase the tally */
                }
                /* reconstruction of flow values */
                idx = idx * 5; /* switch to index field variable */
                U[idx+0] = Uo[0];
                U[idx+1] = Uo[0] * Uo[1];
                U[idx+2] = Uo[0] * Uo[2];
                U[idx+3] = Uo[0] * Uo[3];
                U[idx+4] = Uo[4] / (flow->gamma - 1.0) + 0.5 * Uo[0] * (Uo[1] * Uo[1] + Uo[2] * Uo[2] + Uo[3] * Uo[3]);
            }
        }
    }
    /*
     * Recompute the changed geometry and apply boundary condition.
     */
    ComputeDomainGeometryGCIBM(space, particle, part);
    BoundaryConditionGCIBM(U, space, particle, part, flow);
    return 0;
}
/*
 * Add the pressure force at boundary to the pressure force of corresponding
 * particle. The pressure at boundary will equal to the pressure of the ghost
 * node since zero pressure gradient at wall normal direction is enforced here.
 * A even spaced pressure distribution over the particle surface is assumed 
 * since we only compute the pressure at the boundary point that has a ghost 
 * neighbour. By this approach, the accuracy of pressure integration along 
 * particle surface will increase correspondingly with the increase of mesh 
 * resolution while saving remarkable computation effort.
 */
static int SurfaceForceIntegration(const Real *U, const Space *space, Particle *particle, 
        const Partition *part, const Flow *flow)
{
    int idx = 0; /* linear array index math variable */
    int geoID = 0; /* geometry id */
    Real *ptk = NULL;
    Real distX = 0.0;
    Real distY = 0.0;
    Real distZ = 0.0;
    Real distToCenter = 0.0;
    Real normalX = 0.0;
    Real normalY = 0.0;
    Real normalZ = 0.0;
    Real ds = 0.0;
    const Real pi = acos(-1);
    Real p = 0.0;
    const int offset = space->nodeFlagOffset;
    /* reset some non accumulative information of particles to zero */
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        ptk = particle->headAddress + geoCount * particle->entryN; /* point to storage of current particle */
        ptk[8] = 0; /* force at x direction */
        ptk[9] = 0; /* force at y direction */
        ptk[10] = 0; /* force at z direction */
        ptk[11] = 0; /* ghost node count */
    }
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (offset > space->nodeFlag[idx]) { /* it's not a ghost */
                    continue;
                }
                geoID = space->nodeFlag[idx] - offset; /* extract geometry number from inner ghost node flag */
                ptk = particle->headAddress + geoID * particle->entryN; /* point to storage of current particle */
                distX = space->xMin + (i - space->ng) * space->dx - ptk[0];
                distY = space->yMin + (j - space->ng) * space->dy - ptk[1];
                distZ = space->zMin + (k - space->ng) * space->dz - ptk[2];
                distToCenter = sqrt(distX * distX + distY * distY + distZ * distZ);
                normalX = distX / distToCenter;
                normalY = distY / distToCenter;
                normalZ = distZ / distToCenter;
                idx = idx * 5; /* switch to index field variable */
                p = (flow->gamma - 1.0) * (U[idx+4] - 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx+0]);
                ptk[8] = -p * normalX; /* increase fx by pressure projection on x */
                ptk[9] = -p * normalY; /* increase fy by pressure projection on y */
                ptk[10] = -p * normalZ; /* increase fz by pressure projection on z */
                ptk[11] = ptk[11] + 1; /* count the number of ghost node of current particle */

            }
        }
    }
    /* calibrate the sum of discrete forces into integration */
    if (0 == space->dx * space->dy * space->dz) { /* space dimension collapsed */
        ds = 2 * pi * ptk[3] / ptk[11]; /* circle perimeter */
    } else {
        ds = 4 * pi * ptk[3] * ptk[3] / ptk[11]; /* sphere surface */
    }
    ptk[8] = ptk[8] * ds;
    ptk[9] = ptk[9] * ds;
    ptk[10] = ptk[10] * ds;
    return 0;
}
/* a good practice: end file with a newline */

