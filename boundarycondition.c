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
#include <math.h> /* common mathematical functions */
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Values should be normalized values relative to the reference values.
 */
int BoundaryCondtion(Field *field, const Space *space, const Partition *part, const Flow *flow)
{
    /*
     * Decompose the field variable into each component.
     */
    Real *rho = field->Un + 0 * space->kMax * space->jMax * space->iMax;
    Real *rho_u = field->Un + 1 * space->kMax * space->jMax * space->iMax;
    Real *rho_v = field->Un + 2 * space->kMax * space->jMax * space->iMax;
    Real *rho_w = field->Un + 3 * space->kMax * space->jMax * space->iMax;
    Real *rho_eT = field->Un + 4 * space->kMax * space->jMax * space->iMax;
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int idx = 0; /* calculated index */
    int idxW = 0; /* index at West */
    int idxE = 0; /* index at East */
    int idxS = 0; /* index at South */
    int idxN = 0; /* index at North */
    int idxF = 0; /* index at Front */
    int idxB = 0; /* index at Back */
    int idxWW = 0; /* index at West */
    int idxEE = 0; /* index at East */
    int idxSS = 0; /* index at South */
    int idxNN = 0; /* index at North */
    int idxFF = 0; /* index at Front */
    int idxBB = 0; /* index at Back */
    /*
     * Inlet conditions
     */
    Real rhoInlet = 1.0;
    Real uInlet = 1.0;
    Real vInlet = 0.0;
    Real wInlet = 0.0;
    Real pInlet = 1.0;
    /* Apply conditions to inlet */
    for (k = part->kSub[6]; k < part->kSup[6]; ++k) {
        for (j = part->jSub[6]; j < part->jSup[6]; ++j) {
            for (i = part->iSub[6]; i < part->iSup[6]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                rho[idx] = rhoInlet;
                rho_u[idx] = rhoInlet * uInlet;
                rho_v[idx] = rhoInlet * vInlet;
                rho_w[idx] = rhoInlet * wInlet;
                rho_eT[idx] = pInlet / (flow->gamma - 1) + 
                    0.5 * rhoInlet * (uInlet * uInlet + vInlet * vInlet + wInlet * wInlet);
            }
        }
    }
    /* Extrapolate values to ghost cells by zero gradient condition */
    for (k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (i = part->iSup[0] - 1; i >= part->iSub[0]; --i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxE = (k * space->jMax + j) * space->iMax + i + 1;
                rho[idx] = rho[idxE];
                rho_u[idx] = rho_u[idxE];
                rho_v[idx] = rho_v[idxE];
                rho_w[idx] = rho_w[idxE];
                rho_eT[idx] = rho_eT[idxE];
            }
        }
    }
    /*
     * Outlet conditions
     */
    /* Apply conditions to outlet */
    for (k = part->kSub[7]; k < part->kSup[7]; ++k) {
        for (j = part->jSub[7]; j < part->jSup[7]; ++j) {
            for (i = part->iSub[7]; i < part->iSup[7]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxWW = (k * space->jMax + j) * space->iMax + i - 2;
                rho[idx] = 2 * rho[idxW] - rho[idxWW];
                rho_u[idx] = 2 * rho_u[idxW] - rho_u[idxWW];
                rho_v[idx] = 2 * rho_v[idxW] - rho_v[idxWW];
                rho_w[idx] = 2 * rho_w[idxW] - rho_w[idxWW];
                rho_eT[idx] = 2 * rho_eT[idxW] - rho_eT[idxWW];
            }
        }
    }
    /* Extrapolate values to ghost cells by zero gradient condition */
    for (k = part->kSub[1]; k < part->kSup[1]; ++k) {
        for (j = part->jSub[1]; j < part->jSup[1]; ++j) {
            for (i = part->iSub[1]; i < part->iSup[1]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxWW = (k * space->jMax + j) * space->iMax + i - 2;
                rho[idx] = 2 * rho[idxW] - rho[idxWW];
                rho_u[idx] = 2 * rho_u[idxW] - rho_u[idxWW];
                rho_v[idx] = 2 * rho_v[idxW] - rho_v[idxWW];
                rho_w[idx] = 2 * rho_w[idxW] - rho_w[idxWW];
                rho_eT[idx] = 2 * rho_eT[idxW] - rho_eT[idxWW];
            }
        }
    }
    return 0;
}
/* a good practice: end file with a newline */

