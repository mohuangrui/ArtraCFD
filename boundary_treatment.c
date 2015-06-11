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
static int ZeroGradientFlow(Real *, const int, const int);
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
    Real Uo[DIMUo] = { /* obtain primitive values of current boundary */
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
                        ZeroGradientFlow(U, idx, idxh);
                        break;
                    case 3: /* slip wall, zero-gradient for scalar and tangential component, zero for normal component */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        PrimitiveByConservative(Uoh, idxh, U, flow);
                        Uo[1] = (!normalX) * Uoh[1];
                        Uo[2] = (!normalY) * Uoh[2];
                        Uo[3] = (!normalZ) * Uoh[3];
                        Uo[4] = Uoh[4]; /* zero normal gradient of pressure */
                        if (0 > Uo[5]) { /* adiabatic, dT/dn = 0 */
                            Uo[5] = Uoh[5];
                        } /* otherwise, use specified constant wall temperature, T = Tw */
                        Uo[0] = Uo[4] / (Uo[5] * flow->gasR); /* compute density */
                        ConservativeByPrimitive(U, idx, Uo, flow);
                        break;
                    case 4: /* noslip wall */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        PrimitiveByConservative(Uoh, idxh, U, flow);
                        Uo[1] = 0;
                        Uo[2] = 0;
                        Uo[3] = 0;
                        Uo[4] = Uoh[4]; /* zero normal gradient of pressure */
                        if (0 > Uo[5]) { /* adiabatic, dT/dn = 0 */
                            Uo[5] = Uoh[5];
                        } /* otherwise, use specified constant wall temperature, T = Tw */
                        Uo[0] = Uo[4] / (Uo[5] * flow->gasR); /* compute density */
                        ConservativeByPrimitive(U, idx, Uo, flow);
                        break;
                    case 5: /* primary periodic pair, apply boundary translation */
                        idxh = IndexMath(k - (space->nz - 2) * normalZ, j - (space->ny - 2) * normalY, 
                                i - (space->nx - 2) * normalX, space) * DIMU;
                        ZeroGradientFlow(U, idx, idxh);
                        break;
                    case -5: /* auxiliary periodic pair, apply zero gradient flow */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        ZeroGradientFlow(U, idx, idxh);
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
                            ZeroGradientFlow(U, idx, idxh);
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
                                Uoh[4] / (flow->gamma - 1.0);
                            break;
                    }
                }
            }
        }
    }
    return 0;
}
static int ZeroGradientFlow(Real *U, const int idx, const int idxh)
{
    for (int dim = 0; dim < DIMU; ++dim) {
        U[idx+dim] = U[idxh+dim];
    }
    return 0;
}
/* a good practice: end file with a newline */

