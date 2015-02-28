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
 * Function definitions
 ****************************************************************************/
/*
 * Values should be normalized values relative to the reference values.
 */
int BoundaryCondtion(Real *U, const Space *space, const Particle *particle, 
        const Partition *part, const Flow *flow)
{
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int dim = 0; /* dimension count of vectors */
    int idx = 0; /* linear array index math variable */
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
    const Real rhoInlet = 1.0;
    const Real uInlet = 1.0;
    const Real vInlet = 0.0;
    const Real wInlet = 0.0;
    const Real pInlet = 1.0;
    /* Apply conditions to inlet */
    domainID = 6;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                U[idx+0] = rhoInlet;
                U[idx+1] = rhoInlet * uInlet;
                U[idx+2] = rhoInlet * vInlet;
                U[idx+3] = rhoInlet * wInlet;
                U[idx+4] = pInlet / (flow->gamma - 1) + 
                    0.5 * rhoInlet * (uInlet * uInlet + vInlet * vInlet + wInlet * wInlet);
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 0;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSup[domainID] - 1; i >= part->iSub[domainID]; --i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxE = ((k * space->jMax + j) * space->iMax + i + 1) * 5;
                for (dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = U[idxE+dim];
                }
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
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxW = ((k * space->jMax + j) * space->iMax + i - 1) * 5;
                idxWW = ((k * space->jMax + j) * space->iMax + i - 2) * 5;
                for (dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = 2 * U[idxW+dim] - U[idxWW+dim];
                }
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 1;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxW = ((k * space->jMax + j) * space->iMax + i - 1) * 5;
                idxWW = ((k * space->jMax + j) * space->iMax + i - 2) * 5;
                for (dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = 2 * U[idxW+dim] - U[idxWW+dim];
                }
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
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxN = ((k * space->jMax + j + 1) * space->iMax + i) * 5;
                U[idx+0] = U[idxN+0];
                U[idx+1] = U[idxN+1];
                U[idx+2] = 0;
                U[idx+3] = 0;
                U[idx+4] = U[idxN+4];
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 2;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSup[domainID] - 1; j >= part->jSub[domainID]; --j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxN = ((k * space->jMax + j + 1) * space->iMax + i) * 5;
                idxNN = ((k * space->jMax + j + 2) * space->iMax + i) * 5;
                for (dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = 2 * U[idxN+dim] - U[idxNN+dim];
                }
            }
        }
    }
    /* Apply conditions to north boundary */
    domainID = 9;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxS = ((k * space->jMax + j - 1) * space->iMax + i) * 5;
                U[idx+0] = U[idxS+0];
                U[idx+1] = U[idxS+1];
                U[idx+2] = 0;
                U[idx+3] = 0;
                U[idx+4] = U[idxS+4];
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 3;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxS = ((k * space->jMax + j - 1) * space->iMax + i) * 5;
                idxSS = ((k * space->jMax + j - 2) * space->iMax + i) * 5;
                for (dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = 2 * U[idxS+dim] - U[idxSS+dim];
                }
            }
        }
    }
    /* Apply conditions to front boundary */
    domainID = 10;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxB = (((k + 1) * space->jMax + j) * space->iMax + i) * 5;
                U[idx+0] = U[idxB+0];
                U[idx+1] = U[idxB+1];
                U[idx+2] = 0;
                U[idx+3] = 0;
                U[idx+4] = U[idxB+4];
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 4;
    for (k = part->kSup[domainID] - 1; k >= part->kSub[domainID]; --k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxB = (((k + 1) * space->jMax + j) * space->iMax + i) * 5;
                idxBB = (((k + 2) * space->jMax + j) * space->iMax + i) * 5;
                for (dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = 2 * U[idxB+dim] - U[idxBB+dim];
                }
            }
        }
    }
    /* Apply conditions to back boundary */
    domainID = 11;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxF = (((k - 1) * space->jMax + j) * space->iMax + i) * 5;
                U[idx+0] = U[idxF+0];
                U[idx+1] = U[idxF+1];
                U[idx+2] = 0;
                U[idx+3] = 0;
                U[idx+4] = U[idxF+4];
            }
        }
    }
    /* Extrapolate values to ghost cells */
    domainID = 5;
    for (k = part->kSub[domainID]; k < part->kSup[domainID]; ++k) {
        for (j = part->jSub[domainID]; j < part->jSup[domainID]; ++j) {
            for (i = part->iSub[domainID]; i < part->iSup[domainID]; ++i) {
                idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                idxF = (((k - 1) * space->jMax + j) * space->iMax + i) * 5;
                idxFF = (((k - 2) * space->jMax + j) * space->iMax + i) * 5;
                for (dim = 0; dim < 5; ++dim) {
                    U[idx+dim] = 2 * U[idxF+dim] - U[idxFF+dim];
                }
            }
        }
    }
    /*
     * Boundary condition for interior ghost cells
     */
    BoundaryConditionGCIBM(U, space, particle, part);
    return 0;
}
/* a good practice: end file with a newline */

