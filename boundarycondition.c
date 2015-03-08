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
    BoundaryConditionGCIBM(U, space, particle, part);
    return 0;
}
static int ApplyBoundaryConditions(const int partID, Real *U, const Space *space, 
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
    Real rho_h = 0; 
    Real u_h = 0;
    Real u_hh = 0;
    Real v_h = 0;
    Real v_hh = 0;
    Real w_h = 0;
    Real w_hh = 0;
    Real eT_h = 0;
    Real p_h = 0;
    for (int k = part->kSub[partID]; k < part->kSup[partID]; ++k) {
        for (int j = part->jSub[partID]; j < part->jSup[partID]; ++j) {
            for (int i = part->iSub[partID]; i < part->iSup[partID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                /*
                 * Apply boundary conditions for current node, always remember
                 * that boundary conditions should be based on primitive
                 * variables rather than conservative variables.
                 */
                switch (part->typeBC[partID]) {
                    case 1: /* inflow */
                        U[idx+0] = rho;
                        U[idx+1] = rho * u;
                        U[idx+2] = rho * v;
                        U[idx+3] = rho * w;
                        U[idx+4] = p / (flow->gamma - 1) + 0.5 * rho * (u * u + v * v + w * w);
                        break;
                    case 2: /* outflow */
                        /* Calculate inner neighbour nodes according to normal vector direction. */
                        idxh = (((k - normalZ) * space->jMax + (j - normalY)) * space->iMax + i - normalX) * 5;
                        idxhh = (((k - 2 * normalZ) * space->jMax + (j - 2 * normalY)) * space->iMax + i - 2 * normalX) * 5;
                        rho_h = U[idxh+0];
                        u_h = U[idxh+1] / rho_h;
                        v_h = U[idxh+2] / rho_h;
                        w_h = U[idxh+3] / rho_h;
                        eT_h = U[idxh+4] / rho_h;
                        p_h = (flow->gamma - 1) * rho_h * (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h));
                        u_hh = U[idxhh+1] / U[idxhh+0];
                        v_hh = U[idxhh+2] / U[idxhh+0];
                        w_hh = U[idxhh+3] / U[idxhh+0];
                        U[idx+0] = rho_h;
                        U[idx+1] = rho_h * (2 * u_h - u_hh);
                        U[idx+2] = rho_h * (2 * v_h - v_hh);
                        U[idx+3] = rho_h * (2 * w_h - w_hh);
                        U[idx+4] = p_h / (flow->gamma - 1) + 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx+0];
                        break;
                    case 3: /* slip wall, the method of image, copy scalar, flip vector */
                        idxh = (((k - normalZ) * space->jMax + (j - normalY)) * space->iMax + i - normalX) * 5;
                        rho_h = U[idxh+0];
                        u_h = U[idxh+1] / rho_h;
                        v_h = U[idxh+2] / rho_h;
                        w_h = U[idxh+3] / rho_h;
                        eT_h = U[idxh+4] / rho_h;
                        p_h = (flow->gamma - 1) * rho_h * (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h));
                        U[idx+0] = rho_h;
                        U[idx+1] = (!normalX) * rho_h * u_h;
                        U[idx+2] = (!normalY) * rho_h * v_h;
                        U[idx+3] = (!normalZ) * rho_h * w_h;
                        if (0 > T) { /* adiabatic, dT/dn = 0 */
                            U[idx+4] = p_h / (flow->gamma - 1) + 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx+0];
                        } else { /* constant wall temperature, T = Tw */
                            U[idx+4] = rho_h * flow->cv * T + 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx+0];
                        }
                        break;
                    case 4: /* noslip wall */
                        idxh = (((k - normalZ) * space->jMax + (j - normalY)) * space->iMax + i - normalX) * 5;
                        rho_h = U[idxh+0];
                        u_h = U[idxh+1] / rho_h;
                        v_h = U[idxh+2] / rho_h;
                        w_h = U[idxh+3] / rho_h;
                        eT_h = U[idxh+4] / rho_h;
                        p_h = (flow->gamma - 1) * rho_h * (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h));
                        U[idx+0] = rho_h;
                        U[idx+1] = 0;
                        U[idx+2] = 0;
                        U[idx+3] = 0;
                        if (0 > T) { /* adiabatic, dT/dn = 0 */
                            U[idx+4] = p_h / (flow->gamma - 1);
                        } else { /* constant wall temperature, T = Tw */
                            U[idx+4] = rho_h * flow->cv * T;
                        }
                        break;
                    case 5: /* primary periodic pair, apply boundary translation */
                        idxh = (((k - (space->nz - 2) * normalZ) * space->jMax + (j - (space->ny - 2) * normalY)) * space->iMax + i - (space->nx - 2) * normalX) * 5;
                        for (int dim = 0; dim < 5; ++dim) {
                            U[idx+dim] = U[idxh+dim];
                        }
                        break;
                    case -5: /* auxiliary periodic pair, apply zero gradient flow */
                        idxh = (((k - normalZ) * space->jMax + (j - normalY)) * space->iMax + i - normalX) * 5;
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
                        case 1: /* inflow */
                        case 5: /* primary periodic pair */
                        case -5: /* auxiliary periodic pair */
                            for (int dim = 0; dim < 5; ++dim) {
                                U[idx+dim] = U[idxh+dim];
                            }
                            break;
                        default: /* linear interpolation */
                            rho_h = U[idxh+0];
                            u_h = U[idxh+1] / rho_h;
                            v_h = U[idxh+2] / rho_h;
                            w_h = U[idxh+3] / rho_h;
                            eT_h = U[idxh+4] / rho_h;
                            p_h = (flow->gamma - 1) * rho_h * (eT_h - 0.5 * (u_h * u_h + v_h * v_h + w_h * w_h));
                            u_hh = U[idxhh+1] / U[idxhh+0];
                            v_hh = U[idxhh+2] / U[idxhh+0];
                            w_hh = U[idxhh+3] / U[idxhh+0];
                            U[idx+0] = rho_h;
                            U[idx+1] = rho_h * (2 * u_h - u_hh);
                            U[idx+2] = rho_h * (2 * v_h - v_hh);
                            U[idx+3] = rho_h * (2 * w_h - w_hh);
                            U[idx+4] = p_h / (flow->gamma - 1) + 0.5 * (U[idx+1] * U[idx+1] + U[idx+2] * U[idx+2] + U[idx+3] * U[idx+3]) / U[idx+0];
                            break;
                    }
                }
            }
        }
    }
    return 0;
}
/* a good practice: end file with a newline */

