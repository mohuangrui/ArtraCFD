/****************************************************************************
 * Boundary Condition                                                       *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines the boundary conditions and treatments of the flow.  *
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
static int ApplyBoundaryConditions(const int, Real *, const Space *, 
        const Partition *, const Flow *);
static int ZeroGradientFlow(Real *U, const int idx, const int idxh, const Space *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int BoundaryCondtionsAndTreatments(Real *U, const Space *space, 
        const Particle *particle, const Partition *part, const Flow *flow)
{
    for (int partID = 1; partID < 7; ++partID) {
        ApplyBoundaryConditions(partID, U, space, part, flow);
    }
    /* Boundary conditions and treatments for interior ghost cells */
    BoundaryConditionGCIBM(U, space, particle, part, flow);
    return 0;
}
static int ApplyBoundaryConditions(const int partID, Real *U, const Space *space, 
        const Partition *part, const Flow *flow)
{
    int idx = 0; /* linear array index math variable */
    int idxh = 0; /* index at one node distance */
    int idxhh = 0; /* index at two nodes distance */
    const Real Uo[DIMUo] = { /* obtain primitive values of current boundary */
        part->valueBC[partID][0],
        part->valueBC[partID][1],
        part->valueBC[partID][2],
        part->valueBC[partID][3],
        part->valueBC[partID][4],
        part->valueBC[partID][5]};
    Real Uoh[DIMUo] = {0.0};
    Real Uohh[DIMUo] = {0.0};
    const int normalZ = part->normalZ[partID];
    const int normalY = part->normalY[partID];
    const int normalX = part->normalX[partID];
    for (int k = part->kSub[partID]; k < part->kSup[partID]; ++k) {
        for (int j = part->jSub[partID]; j < part->jSup[partID]; ++j) {
            for (int i = part->iSub[partID]; i < part->iSup[partID]; ++i) {
                idx = IndexMath(k, j, i, space) * DIMU;
                /*
                 * Apply boundary conditions for current node, always remember
                 * that boundary conditions should be based on primitive
                 * variables rather than conservative variables.
                 */
                switch (part->typeBC[partID]) {
                    case 1: /* inflow */
                        ConservativeByPrimitive(U, idx, Uo, flow);
                        break;
                    case 2: /* outflow */
                        /* Calculate inner neighbour nodes according to normal vector direction. */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        ZeroGradientFlow(U, idx, idxh, space);
                        break;
                    case 3: /* slip wall, zero-gradient scalar and tangential component, zero normal component */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        U[idx] = U[idxh];
                        U[idx+1] = (!normalX) * U[idxh+1];
                        U[idx+2] = (!normalY) * U[idxh+2];
                        U[idx+3] = (!normalZ) * U[idxh+3];
                        if (0 > Uo[5]) { /* adiabatic, dT/dn = 0 */
                            U[idx+4] = 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx] + 
                                ComputePressure(idxh, U, flow) / flow->gammaMinusOne;
                        } else { /* constant wall temperature, T = Tw */
                            U[idx+4] = 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx] +
                                U[idx] * Uo[5] * flow->cv;
                        }
                        break;
                    case 4: /* noslip wall */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        U[idx] = U[idxh];
                        U[idx+1] = 0.0;
                        U[idx+2] = 0.0;
                        U[idx+3] = 0.0;
                        if (0 > Uo[5]) { /* adiabatic, dT/dn = 0 */
                            U[idx+4] = ComputePressure(idxh, U, flow) / flow->gammaMinusOne;
                        } else { /* constant wall temperature, T = Tw */
                            U[idx+4] = U[idx] * Uo[5] * flow->cv;
                        }
                        break;
                    case 5: /* primary periodic pair, apply boundary translation */
                        idxh = IndexMath(k - (space->nz - 2) * normalZ, j - (space->ny - 2) * normalY, 
                                i - (space->nx - 2) * normalX, space) * DIMU;
                        ZeroGradientFlow(U, idx, idxh, space);
                        break;
                    case -5: /* auxiliary periodic pair, apply zero gradient flow */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        ZeroGradientFlow(U, idx, idxh, space);
                        break;
                    default:
                        break;
                }
                /*
                 * Extrapolate values for exterior ghost nodes of current node
                 */
                for (int ng = 1; ng <= space->ng; ++ng) { /* process layer by layer */
                    idx = IndexMath(k + ng * normalZ, j + ng * normalY, i + ng * normalX, space) * DIMU;
                    idxh = IndexMath(k + (ng-1) * normalZ, j + (ng-1) * normalY, i + (ng-1) * normalX, space) * DIMU;
                    switch (part->typeBC[partID]) {
                        case 1: /* inflow */
                        case 2: /* outflow */
                        case 5: /* primary periodic pair */
                        case -5: /* auxiliary periodic pair */
                            ZeroGradientFlow(U, idx, idxh, space);
                            break;
                        default: /* linear interpolation will automatically apply the method of image */
                            idxhh = IndexMath(k + (ng-2) * normalZ, j + (ng-2) * normalY, i + (ng-2) * normalX, space) * DIMU;
                            PrimitiveByConservative(Uoh, idxh, U, flow);
                            PrimitiveByConservative(Uohh, idxhh, U, flow);
                            U[idx] = Uoh[0];
                            U[idx+1] = U[idx] * (2.0 * Uoh[1] - Uohh[1]);
                            U[idx+2] = U[idx] * (2.0 * Uoh[2] - Uohh[2]);
                            U[idx+3] = U[idx] * (2.0 * Uoh[3] - Uohh[3]);
                            U[idx+4] = 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx] + 
                                Uoh[4] / flow->gammaMinusOne;
                            break;
                    }
                }
            }
        }
    }
    return 0;
}
static int ZeroGradientFlow(Real *U, const int idx, const int idxh, const Space *space)
{
    /* no loop to avoid loop averload */
    U[idx] = U[idxh];
    U[idx+1] = U[idxh+1];
    U[idx+2] = U[idxh+2];
    U[idx+3] = U[idxh+3];
    U[idx+4] = U[idxh+4];
    return 0;
}
/* a good practice: end file with a newline */

