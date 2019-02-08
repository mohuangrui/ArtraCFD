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
static void ApplyBoundaryCondition(const int, const int, int [restrict][LIMIT],
        const int, Space *, const Model *);
static void EnforceZeroGradient(const Real [restrict], Real [restrict]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void TreatBoundary(const int tn, Space *space, const Model *model)
{
    /*
     * Internal boundary treatment
     * Should be performed first to ensure stencils for diffusive flux
     * discretization are treated correctly, especially for collapsed
     * dimensions.
     */
    TreatImmersedBoundary(tn, space, model);
    /*
     * External boundary treatment
     * When no mixed derivatives are discretized, only cross-type stencils
     * are needed. Then, the corner ghost nodes do not need treatment.
     * However, corner ghost nodes are necessary for the computation of
     * mixed derivatives as in viscous fluxes or for transfer operators
     * within multigrid. To treat the entire ghost region, the boundary
     * treatment should be performed one box layer by one box layer from
     * inside to outside.
     */
    const Partition *const part = &(space->part);
    const IntVec ng = {part->ng[X], part->ng[Y], part->ng[Z]};
    const int R = MaxInt(ng[X], MaxInt(ng[Y], ng[Z]));
    int box[DIMS][LIMIT] = {{0}}; /* range box of numerical boundary */
    for (int r = 0; r <= R; ++r) { /* process layer by layer */
        for (int p = PWB; p <= PBB; ++p) {
            const IntVec N = {part->N[p][X], part->N[p][Y], part->N[p][Z]};
            for (int s = 0; s < DIMS; ++s) { /* compute range box of each layer */
                box[s][MIN] = part->ns[p][s][MIN] + MinInt(r, ng[s]) * (N[s] - !N[s]);
                box[s][MAX] = part->ns[p][s][MAX] + MinInt(r, ng[s]) * (N[s] + !N[s]) - (ng[s] < r) * (!!N[s]);
            }
            if ((box[X][MIN] >= box[X][MAX]) || (box[Y][MIN] >= box[Y][MAX]) || (box[Z][MIN] >= box[Z][MAX])) {
                continue;
            }
            ApplyBoundaryCondition(p, r, box, tn, space, model);
        }
    }
    return;
}
static void ApplyBoundaryCondition(const int p, const int r, int box[restrict][LIMIT],
        const int tn, Space *space, const Model *model)
{
    const Partition *const part = &(space->part);
    Node *const node = space->node;
    const Real zero = 0.0;
    const Real UoGiven[DIMUo] = { /* specified primitive values of current boundary */
        part->varBC[p][0],
        part->varBC[p][1],
        part->varBC[p][2],
        part->varBC[p][3],
        part->varBC[p][4],
        part->varBC[p][5]};
    const IntVec N = {part->N[p][X], part->N[p][Y], part->N[p][Z]};
    const IntVec LN = {part->m[X] * N[X], part->m[Y] * N[Y], part->m[Z] * N[Z]};
    Real *restrict UG = NULL;
    Real *restrict UI = NULL;
    Real *restrict UO = NULL;
    Real *restrict Uh = NULL;
    int idxG = 0; /* index at ghost node */
    int idxI = 0; /* index at image node */
    int idxO = 0; /* index at boundary point */
    int idxh = 0; /* index at neighbouring point */
    Real UoG[DIMUo] = {zero};
    Real UoI[DIMUo] = {zero};
    Real UoO[DIMUo] = {zero};
    Real Uoh[DIMUo] = {zero};
    for (int k = box[Z][MIN]; k < box[Z][MAX]; ++k) {
        for (int j = box[Y][MIN]; j < box[Y][MAX]; ++j) {
            for (int i = box[X][MIN]; i < box[X][MAX]; ++i) {
                /*
                 * Apply boundary conditions for current node, always remember
                 * that boundary conditions should be based on primitive
                 * variables rather than conservative variables.
                 */
                if (0 != r) { /* treat ghost layers */
                    idxG = IndexNode(k, j, i, part->n[Y], part->n[X]);
                    UG = node[idxG].U[tn];
                    switch (part->typeBC[p]) {
                        case SLIPWALL:
                            /* fall through */
                        case NOSLIPWALL:
                            idxO = IndexNode(k - r*N[Z], j - r*N[Y], i - r*N[X], part->n[Y], part->n[X]);
                            UO = node[idxO].U[tn];
                            MapPrimitive(model->gamma, model->gasR, UO, UoO);
                            idxI = IndexNode(k - 2*r*N[Z], j - 2*r*N[Y], i - 2*r*N[X], part->n[Y], part->n[X]);
                            UI = node[idxI].U[tn];
                            MapPrimitive(model->gamma, model->gasR, UI, UoI);
                            DoMethodOfImage(UoI, UoO, UoG);
                            UoG[0] = UoG[4] / (UoG[5] * model->gasR); /* compute density */
                            MapConservative(model->gamma, UoG, UG);
                            break;
                        case PERIODIC:
                            idxh = IndexNode(k - LN[Z], j - LN[Y], i - LN[X], part->n[Y], part->n[X]);
                            Uh = node[idxh].U[tn];
                            EnforceZeroGradient(Uh, UG);
                            break;
                        default:
                            idxh = IndexNode(k - N[Z], j - N[Y], i - N[X], part->n[Y], part->n[X]);
                            Uh = node[idxh].U[tn];
                            EnforceZeroGradient(Uh, UG);
                            break;
                    }
                    continue;
                }
                idxO = IndexNode(k, j, i, part->n[Y], part->n[X]);
                UO = node[idxO].U[tn];
                switch (part->typeBC[p]) { /* treat physical boundary */
                    case INFLOW:
                        MapConservative(model->gamma, UoGiven, UO);
                        break;
                    case OUTFLOW:
                        /* Calculate inner neighbour nodes according to normal vector direction. */
                        idxh = IndexNode(k - N[Z], j - N[Y], i - N[X], part->n[Y], part->n[X]);
                        Uh = node[idxh].U[tn];
                        EnforceZeroGradient(Uh, UO);
                        break;
                    case SLIPWALL: /* zero-gradient for scalar and tangential component, zero for normal component */
                        idxh = IndexNode(k - N[Z], j - N[Y], i - N[X], part->n[Y], part->n[X]);
                        Uh = node[idxh].U[tn];
                        MapPrimitive(model->gamma, model->gasR, Uh, Uoh);
                        UoO[1] = (!N[X]) * Uoh[1];
                        UoO[2] = (!N[Y]) * Uoh[2];
                        UoO[3] = (!N[Z]) * Uoh[3];
                        UoO[4] = Uoh[4]; /* zero normal gradient of pressure */
                        if (zero > UoGiven[5]) { /* adiabatic, dT/dn = 0 */
                            UoO[5] = Uoh[5];
                        } else { /* otherwise, use specified constant wall temperature, T = Tw */
                            UoO[5] = UoGiven[5];
                        }
                        UoO[0] = UoO[4] / (UoO[5] * model->gasR); /* compute density */
                        MapConservative(model->gamma, UoO, UO);
                        break;
                    case NOSLIPWALL:
                        idxh = IndexNode(k - N[Z], j - N[Y], i - N[X], part->n[Y], part->n[X]);
                        Uh = node[idxh].U[tn];
                        MapPrimitive(model->gamma, model->gasR, Uh, Uoh);
                        UoO[1] = zero;
                        UoO[2] = zero;
                        UoO[3] = zero;
                        UoO[4] = Uoh[4]; /* zero normal gradient of pressure */
                        if (zero > UoGiven[5]) { /* adiabatic, dT/dn = 0 */
                            UoO[5] = Uoh[5];
                        } else { /* otherwise, use specified constant wall temperature, T = Tw */
                            UoO[5] = UoGiven[5];
                        }
                        UoO[0] = UoO[4] / (UoO[5] * model->gasR); /* compute density */
                        MapConservative(model->gamma, UoO, UO);
                        break;
                    case PERIODIC:
                        /* no treatment needed since the boundary participates normal computation */
                        break;
                    default:
                        break;
                }
            }
        }
    }
    return;
}
static void EnforceZeroGradient(const Real Uh[restrict], Real U[restrict])
{
    for (int n = 0; n < DIMU; ++n) {
        U[n] = Uh[n];
    }
    return;
}
/* a good practice: end file with a newline */

