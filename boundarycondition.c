/****************************************************************************
 * Boundary Condition                                                       *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines the boundary conditions of the flow.                 *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "boundarycondition.h"
#include <stdio.h> /* standard library for input and output */
#include "gcibm.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ApplyBoundaryCondition(const int, Real *, const Space *, 
        const Partition *, const Flow *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int BoundaryCondtion(Real *U, const Space *space, const Particle *particle, 
        const Partition *part, const Flow *flow)
{
    for (int partID = 1; partID < 7; ++partID) {
        ApplyBoundaryCondition(partID, U, space, part, flow);
    }
    /*
     * Boundary condition for interior ghost cells
     */
    BoundaryConditionGCIBM(U, space, particle, part);
    return 0;
}
static int ApplyBoundaryCondition(const int partID, Real *U, const Space *space, 
        const Partition *part, const Flow *flow)
{
    int idx = 0; /* linear array index math variable */
    int idxh = 0; /* index at one node distance */
    int idxhh = 0; /* index at two nodes distance */
    const Real rho = part->valueBC[partID][0];
    const Real u = part->valueBC[partID][1];
    const Real v = part->valueBC[partID][2];
    const Real w = part->valueBC[partID][3];
    const Real p = part->valueBC[partID][4];
    const Real T = part->valueBC[partID][5];
    const int normalZ = part->normalZ[partID];
    const int normalY = part->normalY[partID];
    const int normalX = part->normalX[partID];
    for (int k = part->kSub[partID]; k < part->kSup[partID]; ++k) {
        for (int j = part->jSub[partID]; j < part->jSup[partID]; ++j) {
            for (int i = part->iSub[partID]; i < part->iSup[partID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                /*
                 * Calculate inner neighbour nodes according to normal vector direction.
                 */
                idxh = (((k - normalZ) * space->jMax + (j - normalY)) * space->iMax + i - normalX) * 5;
                idxhh = (((k - 2 * normalZ) * space->jMax + (j - 2 * normalY)) * space->iMax + i - 2 * normalX) * 5;
                /*
                 * apply boundary condition for current node
                 */
                switch (part->typeBC[partID]) {
                    case 1: /* inlet */
                        U[idx+0] = rho;
                        U[idx+1] = rho * u;
                        U[idx+2] = rho * v;
                        U[idx+3] = rho * w;
                        U[idx+4] = p / (flow->gamma - 1) + 0.5 * rho * (u * u + v * v + w * w);
                        break;
                    case 2: /* outflow */
                        for (int dim = 0; dim < 5; ++dim) {
                            U[idx+dim] = 2 * U[idxh+dim] - U[idxhh+dim];
                        }
                        break;
                    case 3: /* slip wall */
                        U[idx+0] = U[idxh+0];
                        U[idx+1] = (!normalX) * U[idxh+1];
                        U[idx+2] = (!normalY) * U[idxh+2];
                        U[idx+3] = (!normalZ) * U[idxh+3];
                        if (T < 0) { /* adiabatic */
                            U[idx+4] = U[idxh+4];
                        } else {
                            U[idx+4] = U[idx+0] * flow->cv * T + 
                                0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx+0];
                        }
                        break;
                    case 4: /* nonslip wall */
                        U[idx+0] = U[idxh+0];
                        U[idx+1] = 0;
                        U[idx+2] = 0;
                        U[idx+3] = 0;
                        if (T < 0) { /* adiabatic */
                            U[idx+4] = U[idxh+4];
                        } else {
                            U[idx+4] = U[idx+0] * flow->cv * T;
                        }
                        break;
                    case 5: /* primary periodic pair, apply boundary translation */
                        idxh = (((k - (space->nz - 2) * normalZ) * space->jMax + (j - (space->ny - 2) * normalY)) * space->iMax + i - (space->nx - 2) * normalX) * 5;
                        for (int dim = 0; dim < 5; ++dim) {
                            U[idx+dim] = U[idxh+dim];
                        }
                        break;
                    case -5: /* auxiliary periodic pair, apply zero gradient flow */
                        for (int dim = 0; dim < 5; ++dim) {
                            U[idx+dim] = U[idxh+dim];
                        }
                        break;
                    default:
                        break;
                }
                /*
                 * Extrapolate values for exterior ghost nodes of current node
                 */
                for (int ng = 1; ng <= space->ng; ++ng) {
                    idx = (((k + ng * normalZ) * space->jMax + (j + ng * normalY)) * space->iMax + i + ng * normalX) * 5;
                    idxh = (((k + (ng-1) * normalZ) * space->jMax + (j + (ng-1) * normalY)) * space->iMax + i + (ng-1) * normalX) * 5;
                    idxhh = (((k + (ng-2) * normalZ) * space->jMax + (j + (ng-2) * normalY)) * space->iMax + i + (ng-2) * normalX) * 5;
                    switch (part->typeBC[partID]) {
                        case 1: /* inlet */
                        case 5: /* primary periodic pair */
                        case -5: /* auxiliary periodic pair */
                            for (int dim = 0; dim < 5; ++dim) {
                                U[idx+dim] = U[idxh+dim];
                            }
                            break;
                        default: /* linear interpolation */
                            for (int dim = 0; dim < 5; ++dim) {
                                U[idx+dim] = 2 * U[idxh+dim] - U[idxhh+dim];
                            }
                            break;
                    }
                }
            }
        }
    }
    return 0;
}
/* a good practice: end file with a newline */

