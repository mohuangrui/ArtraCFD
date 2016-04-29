/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>                        *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "boundary_treatment.h"
#include <stdio.h> /* standard library for input and output */
#include "immersed_boundary.h"
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static void ApplyBoundaryConditions(const int, const int, Space *, const Model *);
static int ZeroGradient(Real *, const int, const int);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void BoundaryCondtionsAndTreatments(const int tn, Space *space, const Model *model)
{
    for (int p = PWB; p < PWG; ++p) {
        ApplyBoundaryConditions(p, tn, space, model);
    }
    ImmersedBoundaryTreatment(tn, space, model);
    return;
}
static void ApplyBoundaryConditions(const int p, const int tn, 
        Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Node *node = space->node;
    Real *restrict UG = NULL;
    Real *restrict UI = NULL;
    Real *restrict UO = NULL;
    int idxG = 0; /* index at ghost node */
    int idxI = 0; /* index at image node */
    int idxO = 0; /* index at boundary point */
    Real UoG[DIMUo] = {0.0};
    Real UoI[DIMUo] = {0.0};
    Real UoO[DIMUo] = {0.0};
    const Real UoGiven[DIMUo] = { /* specified primitive values of current boundary */
        part->valueBC[p][0],
        part->valueBC[p][1],
        part->valueBC[p][2],
        part->valueBC[p][3],
        part->valueBC[p][4],
        part->valueBC[p][5]};
    const IntVector nl = {part->normal[p][X], part->normal[p][Y], part->normal[p][Z]};
    for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
        for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
            for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                idxO = IndexNode(k, j, i, part->n[Y], part->n[X]);
                U = node[idx].U[TO];
                /*
                 * Apply boundary conditions for current node, always remember
                 * that boundary conditions should be based on primitive
                 * variables rather than conservative variables.
                 */
                switch (part->typeBC[p]) {
                    case INFLOW:
                        ConservativeByPrimitive(U, idxBC, UoGiven, model);
                        break;
                    case 2: /* outflow */
                        /* Calculate inner neighbour nodes according to normal vector direction. */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        ZeroGradient(U, idxBC, idxh);
                        break;
                    case 3: /* slip wall, zero-gradient for scalar and tangential component, zero for normal component */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        PrimitiveByConservative(Uoh, idxh, U, model);
                        UoBC[1] = (!normalX) * Uoh[1];
                        UoBC[2] = (!normalY) * Uoh[2];
                        UoBC[3] = (!normalZ) * Uoh[3];
                        UoBC[4] = Uoh[4]; /* zero normal gradient of pressure */
                        if (0 > UoGiven[5]) { /* adiabatic, dT/dn = 0 */
                            UoBC[5] = Uoh[5];
                        } else { /* otherwise, use specified constant wall temperature, T = Tw */
                            UoBC[5] = UoGiven[5];
                        }
                        UoBC[0] = UoBC[4] / (UoBC[5] * model->gasR); /* compute density */
                        ConservativeByPrimitive(U, idxBC, UoBC, model);
                        break;
                    case 4: /* noslip wall */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        PrimitiveByConservative(Uoh, idxh, U, model);
                        UoBC[1] = 0;
                        UoBC[2] = 0;
                        UoBC[3] = 0;
                        UoBC[4] = Uoh[4]; /* zero normal gradient of pressure */
                        if (0 > UoGiven[5]) { /* adiabatic, dT/dn = 0 */
                            UoBC[5] = Uoh[5];
                        } else { /* otherwise, use specified constant wall temperature, T = Tw */
                            UoBC[5] = UoGiven[5];
                        }
                        UoBC[0] = UoBC[4] / (UoBC[5] * model->gasR); /* compute density */
                        ConservativeByPrimitive(U, idxBC, UoBC, model);
                        break;
                    case 5: /* primary periodic pair, apply boundary translation */
                        idxh = IndexMath(k - (space->nz - 2) * normalZ, j - (space->ny - 2) * normalY, 
                                i - (space->nx - 2) * normalX, space) * DIMU;
                        ZeroGradient(U, idxBC, idxh);
                        break;
                    case -5: /* auxiliary periodic pair, apply zero gradient */
                        idxh = IndexMath(k - normalZ, j - normalY, i - normalX, space) * DIMU;
                        ZeroGradient(U, idxBC, idxh);
                        break;
                    default:
                        break;
                }
                /*
                 * Reconstruct values for exterior ghost nodes of current node
                 */
                for (int ng = 1; ng <= space->ng; ++ng) { /* process layer by layer */
                    idxGhost = IndexMath(k + ng * normalZ, j + ng * normalY, i + ng * normalX, space) * DIMU;
                    switch (part->typeBC[p]) {
                        case 3: /* slip wall */
                        case 4: /* noslip wall */
                            /* 
                             * Apply the method of image.
                             *  -- reflecting vectors over wall for both slip and noslip, stationary and
                             *     moving conditions is unified by linear interpolation.
                             *  -- scalars are symmetrically reflected between image and ghost.
                             */
                            idxImage = IndexMath(k - ng * normalZ, j - ng * normalY, i - ng * normalX, space) * DIMU;
                            PrimitiveByConservative(UoBC, idxBC, U, model);
                            PrimitiveByConservative(UoImage, idxImage, U, model);
                            UoGhost[1] = 2.0 * UoBC[1] - UoImage[1];
                            UoGhost[2] = 2.0 * UoBC[2] - UoImage[2];
                            UoGhost[3] = 2.0 * UoBC[3] - UoImage[3];
                            UoGhost[4] = UoImage[4];
                            UoGhost[5] = UoImage[5];
                            UoGhost[0] = UoGhost[4] / (UoGhost[5] * model->gasR); /* compute density */
                            ConservativeByPrimitive(U, idxGhost, UoGhost, model);
                            break;
                        default:
                            ZeroGradient(U, idxGhost, idxBC);
                            break;
                    }
                }
            }
        }
    }
    return 0;
}
static int ZeroGradient(Real *U, const int idx, const int idxh)
{
    for (int dim = 0; dim < DIMU; ++dim) {
        U[idx+dim] = U[idxh+dim];
    }
    return 0;
}
/* a good practice: end file with a newline */

