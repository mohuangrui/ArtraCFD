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
    Real *Un[5] = {
        field->Un + 0 * space->nMax,
        field->Un + 1 * space->nMax,
        field->Un + 2 * space->nMax,
        field->Un + 3 * space->nMax,
        field->Un + 4 * space->nMax};
    /*
     * Indices
     */
    int k = 0; /* loop count */
    int j = 0; /* loop count */
    int i = 0; /* loop count */
    int dim = 0; /* dimension count of vectors */
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
                Un[0][idx] = rhoInlet;
                Un[1][idx] = rhoInlet * uInlet;
                Un[2][idx] = rhoInlet * vInlet;
                Un[3][idx] = rhoInlet * wInlet;
                Un[4][idx] = pInlet / (flow->gamma - 1) + 
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
                for (dim = 0; dim < 5; ++dim) {
                    Un[dim][idx] = Un[dim][idxE];
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
                idx = (k * space->jMax + j) * space->iMax + i;
                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxWW = (k * space->jMax + j) * space->iMax + i - 2;
                for (dim = 0; dim < 5; ++dim) {
                    Un[dim][idx] = 2 * Un[dim][idxW] - Un[dim][idxWW];
                }
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
                for (dim = 0; dim < 5; ++dim) {
                    Un[dim][idx] = 2 * Un[dim][idxW] - Un[dim][idxWW];
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
                idx = (k * space->jMax + j) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                Un[0][idx] = Un[0][idxN];
                Un[1][idx] = Un[1][idxN];
                Un[2][idx] = 0;
                Un[3][idx] = 0;
                Un[4][idx] = Un[4][idxN];
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
                for (dim = 0; dim < 5; ++dim) {
                    Un[dim][idx] = 2 * Un[dim][idxN] - Un[dim][idxNN];
                }
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
                Un[0][idx] = Un[0][idxS];
                Un[1][idx] = Un[1][idxS];
                Un[2][idx] = 0;
                Un[3][idx] = 0;
                Un[4][idx] = Un[4][idxS];
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
                for (dim = 0; dim < 5; ++dim) {
                    Un[dim][idx] = 2 * Un[dim][idxS] - Un[dim][idxSS];
                }
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
                Un[0][idx] = Un[0][idxB];
                Un[1][idx] = Un[1][idxB];
                Un[2][idx] = 0;
                Un[3][idx] = 0;
                Un[4][idx] = Un[4][idxB];
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
                for (dim = 0; dim < 5; ++dim) {
                    Un[dim][idx] = 2 * Un[dim][idxB] - Un[dim][idxBB];
                }
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
                Un[0][idx] = Un[0][idxF];
                Un[1][idx] = Un[1][idxF];
                Un[2][idx] = 0;
                Un[3][idx] = 0;
                Un[4][idx] = Un[4][idxF];
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
                for (dim = 0; dim < 5; ++dim) {
                    Un[dim][idx] = 2 * Un[dim][idxF] - Un[dim][idxFF];
                }
            }
        }
    }
    /*
     * Wall boundary condition for interior ghost cells
     */
    for (k = part->kSub[12]; k < part->kSup[12]; ++k) {
        for (j = part->jSub[12]; j < part->jSup[12]; ++j) {
            for (i = part->iSub[12]; i < part->iSup[12]; ++i) {
                idx = (k * space->jMax + j) * space->iMax + i;
                if (space->ghostFlag[idx] != 1) { /* it's not a ghost */
                    continue;
                }
                idxW = (k * space->jMax + j) * space->iMax + i - 1;
                idxE = (k * space->jMax + j) * space->iMax + i + 1;
                idxS = (k * space->jMax + j - 1) * space->iMax + i;
                idxN = (k * space->jMax + j + 1) * space->iMax + i;
                idxF = ((k - 1) * space->jMax + j) * space->iMax + i;
                idxB = ((k + 1) * space->jMax + j) * space->iMax + i;
            }
        }
    }
    return 0;
}
/* a good practice: end file with a newline */

