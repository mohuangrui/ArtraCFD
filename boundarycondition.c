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
    Real *rho = field->Un + 0 * space->nMax;
    Real *rho_u = field->Un + 1 * space->nMax;
    Real *rho_v = field->Un + 2 * space->nMax;
    Real *rho_w = field->Un + 3 * space->nMax;
    Real *rho_eT = field->Un + 4 * space->nMax;
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
    int idxSS = 0; /* index at South */
    int idxNN = 0; /* index at North */
    int idxFF = 0; /* index at Front */
    int idxBB = 0; /* index at Back */
    int domainID = 0; /* ID of domain */
    /*
     * Inlet conditions
     */
    Real rhoInlet = 1.0;
    Real uInlet = 1.0;
    Real vInlet = 0.0;
    Real wInlet = 0.0;
    Real pInlet = 1.0;
    /* Apply conditions to inlet */
    domainID = 6;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
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
    /* Extrapolate values to ghost cells */
    domainID = 0;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSup[domainID] - 1; i >= part->iSub[domainID]; --i) {
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
    domainID = 7;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
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
    /* Extrapolate values to ghost cells */
    domainID = 1;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
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
    /*
     * Wall conditions
     */
    /* Apply conditions to south boundary */
    domainID = 8;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                rho[idx] = rho[idxN];
                rho_u[idx] = rho_u[idxN];
                rho_v[idx] = 0;
                rho_w[idx] = 0;
                rho_eT[idx] = rho_eT[idxN];
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 2;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSup[domainID] - 1; j >= part->jSub[domainID]; --j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                idxNN = (k * space->jMax + j + 2) * space->iMax + i;
                rho[idx] = 2 * rho[idxN] - rho[idxNN];
                rho_u[idx] = 2 * rho_u[idxN] - rho_u[idxNN];
                rho_v[idx] = 2 * rho_v[idxN] - rho_v[idxNN];
                rho_w[idx] = 2 * rho_w[idxN] - rho_w[idxNN];
                rho_eT[idx] = 2 * rho_eT[idxN] - rho_eT[idxNN];
            }
        }
    }
    /* Apply conditions to north boundary */
    domainID = 9;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxS = (k * space->jMax + j - 1) * space->iMax + i;
                rho[idx] = rho[idxS];
                rho_u[idx] = rho_u[idxS];
                rho_v[idx] = 0;
                rho_w[idx] = 0;
                rho_eT[idx] = rho_eT[idxS];
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 3;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxS = (k * space->jMax + j - 1) * space->iMax + i;
                idxSS = (k * space->jMax + j - 2) * space->iMax + i;
                rho[idx] = 2 * rho[idxS] - rho[idxSS];
                rho_u[idx] = 2 * rho_u[idxS] - rho_u[idxSS];
                rho_v[idx] = 2 * rho_v[idxS] - rho_v[idxSS];
                rho_w[idx] = 2 * rho_w[idxS] - rho_w[idxSS];
                rho_eT[idx] = 2 * rho_eT[idxS] - rho_eT[idxSS];
            }
        }
    }
    /* Apply conditions to front boundary */
    domainID = 10;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxB = ((k + 1) * space->jMax + j) * space->iMax + i;
                rho[idx] = rho[idxB];
                rho_u[idx] = rho_u[idxB];
                rho_v[idx] = 0;
                rho_w[idx] = 0;
                rho_eT[idx] = rho_eT[idxB];
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 4;
    for (k = part->kSup[domainID] - 1; k >= part->kSub[domainID]; --k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxB = ((k + 1) * space->jMax + j) * space->iMax + i;
                idxBB = ((k + 2) * space->jMax + j) * space->iMax + i;
                rho[idx] = 2 * rho[idxB] - rho[idxBB];
                rho_u[idx] = 2 * rho_u[idxB] - rho_u[idxBB];
                rho_v[idx] = 2 * rho_v[idxB] - rho_v[idxBB];
                rho_w[idx] = 2 * rho_w[idxB] - rho_w[idxBB];
                rho_eT[idx] = 2 * rho_eT[idxB] - rho_eT[idxBB];
            }
        }
    }
    /* Apply conditions to back boundary */
    domainID = 11;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxF = ((k - 1) * space->jMax + j) * space->iMax + i;
                rho[idx] = rho[idxF];
                rho_u[idx] = rho_u[idxF];
                rho_v[idx] = 0;
                rho_w[idx] = 0;
                rho_eT[idx] = rho_eT[idxF];
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 5;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                idxF = ((k - 1) * space->jMax + j) * space->iMax + i;
                idxFF = ((k - 2) * space->jMax + j) * space->iMax + i;
                rho[idx] = 2 * rho[idxF] - rho[idxFF];
                rho_u[idx] = 2 * rho_u[idxF] - rho_u[idxFF];
                rho_v[idx] = 2 * rho_v[idxF] - rho_v[idxFF];
                rho_w[idx] = 2 * rho_w[idxF] - rho_w[idxFF];
                rho_eT[idx] = 2 * rho_eT[idxF] - rho_eT[idxFF];
            }
        }
    }
    return 0;
}
/* a good practice: end file with a newline */

