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
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Add the pressure force at boundary to the pressure force of 
 * corresponding particle. The pressure at boundary will equal
 * to the pressure of the ghost node since zero pressure 
 * gradient at wall normal direction is enforced here.
 * A even spaced pressure distribution over the particle 
 * surface is assumed since we only compute the pressure at
 * the boundary point that has a ghost neighbour. By this 
 * approach, the accuracy of pressure integration along particle
 * surface will increase correspondingly with the increase of 
 * mesh resolution while saving remarkable computation effort.
 */
int SurfaceForceIntegration(const Real *U, const Space *space, Particle *particle, 
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
}
/* a good practice: end file with a newline */

