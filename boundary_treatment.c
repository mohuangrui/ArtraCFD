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
static void ZeroGradient(const Real [restrict], Real [restrict]);
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
static void ApplyBoundaryConditions(const int p, const int tn, Space *space, const Model *model)
{
    const Partition *restrict part = &(space->part);
    Node *const node = space->node;
    Real *restrict UG = NULL;
    Real *restrict UI = NULL;
    Real *restrict UO = NULL;
    Real *restrict Uh = NULL;
    int idxG = 0; /* index at ghost node */
    int idxI = 0; /* index at image node */
    int idxO = 0; /* index at boundary point */
    int idxh = 0; /* index at neighbouring point */
    Real UoG[DIMUo] = {0.0};
    Real UoI[DIMUo] = {0.0};
    Real UoO[DIMUo] = {0.0};
    Real Uoh[DIMUo] = {0.0};
    const Real UoGiven[DIMUo] = { /* specified primitive values of current boundary */
        part->valueBC[p][0],
        part->valueBC[p][1],
        part->valueBC[p][2],
        part->valueBC[p][3],
        part->valueBC[p][4],
        part->valueBC[p][5]};
    const IntVec N = {part->N[p][X], part->N[p][Y], part->N[p][Z]};
    const IntVec LN = {(part->m[X] - 1) * N[X], (part->m[Y] - 1) * N[Y], (part->m[Z] - 1) * N[Z]};
    for (int k = part->ns[p][Z][MIN]; k < part->ns[p][Z][MAX]; ++k) {
        for (int j = part->ns[p][Y][MIN]; j < part->ns[p][Y][MAX]; ++j) {
            for (int i = part->ns[p][X][MIN]; i < part->ns[p][X][MAX]; ++i) {
                idxO = IndexNode(k, j, i, part->n[Y], part->n[X]);
                UO = node[idxO].U[tn];
                /*
                 * Apply boundary conditions for current node, always remember
                 * that boundary conditions should be based on primitive
                 * variables rather than conservative variables.
                 */
                switch (part->typeBC[p]) {
                    case INFLOW:
                        ConservativeByPrimitive(model->gamma, UoGiven, UO);
                        break;
                    case OUTFLOW:
                        /* Calculate inner neighbour nodes according to normal vector direction. */
                        idxh = IndexNode(k - N[Z], j - N[Y], i - N[X], part->n[Y], part->n[X]);
                        Uh = node[idxh].U[tn];
                        ZeroGradient(Uh, UO);
                        break;
                    case SLIPWALL: /* zero-gradient for scalar and tangential component, zero for normal component */
                        idxh = IndexNode(k - N[Z], j - N[Y], i - N[X], part->n[Y], part->n[X]);
                        Uh = node[idxh].U[tn];
                        PrimitiveByConservative(model->gamma, model->gasR, Uh, Uoh);
                        UoO[1] = (!N[X]) * Uoh[1];
                        UoO[2] = (!N[Y]) * Uoh[2];
                        UoO[3] = (!N[Z]) * Uoh[3];
                        UoO[4] = Uoh[4]; /* zero normal gradient of pressure */
                        if (0.0 > UoGiven[5]) { /* adiabatic, dT/dn = 0 */
                            UoO[5] = Uoh[5];
                        } else { /* otherwise, use specified constant wall temperature, T = Tw */
                            UoO[5] = UoGiven[5];
                        }
                        UoO[0] = UoO[4] / (UoO[5] * model->gasR); /* compute density */
                        ConservativeByPrimitive(model->gamma, UoO, UO);
                        break;
                    case NOSLIPWALL:
                        idxh = IndexNode(k - N[Z], j - N[Y], i - N[X], part->n[Y], part->n[X]);
                        Uh = node[idxh].U[tn];
                        PrimitiveByConservative(model->gamma, model->gasR, Uh, Uoh);
                        UoO[1] = 0.0;
                        UoO[2] = 0.0;
                        UoO[3] = 0.0;
                        UoO[4] = Uoh[4]; /* zero normal gradient of pressure */
                        if (0.0 > UoGiven[5]) { /* adiabatic, dT/dn = 0 */
                            UoO[5] = Uoh[5];
                        } else { /* otherwise, use specified constant wall temperature, T = Tw */
                            UoO[5] = UoGiven[5];
                        }
                        UoO[0] = UoO[4] / (UoO[5] * model->gasR); /* compute density */
                        ConservativeByPrimitive(model->gamma, UoO, UO);
                        break;
                    case PERIODIC:
                        idxh = IndexNode(k - LN[Z], j - LN[Y], i - LN[X], part->n[Y], part->n[X]);
                        Uh = node[idxh].U[tn];
                        ZeroGradient(Uh, UO);
                        break;
                    default:
                        break;
                }
                /*
                 * Reconstruct values for exterior ghost nodes of current node
                 */
                for (int ng = 1; ng <= part->ng; ++ng) { /* process layer by layer */
                    idxG = IndexNode(k + ng * N[Z], j + ng * N[Y], i + ng * N[X], part->n[Y], part->n[X]);
                    UG = node[idxG].U[tn];
                    switch (part->typeBC[p]) {
                        case SLIPWALL:
                        case NOSLIPWALL:
                            idxI = IndexNode(k - ng * N[Z], j - ng * N[Y], i - ng * N[X], part->n[Y], part->n[X]);
                            UI = node[idxI].U[tn];
                            PrimitiveByConservative(model->gamma, model->gasR, UI, UoI);
                            MethodOfImage(UoI, UoO, UoG);
                            UoG[0] = UoG[4] / (UoG[5] * model->gasR); /* compute density */
                            ConservativeByPrimitive(model->gamma, UoG, UG);
                            break;
                        case PERIODIC:
                            idxh = IndexNode(k + ng * N[Z] - LN[Z], j + ng * N[Y] - LN[Y], i + ng * N[X] - LN[X], 
                                    part->n[Y], part->n[X]);
                            Uh = node[idxh].U[tn];
                            ZeroGradient(Uh, UG);
                            break;
                        default:
                            ZeroGradient(UO, UG);
                            break;
                    }
                }
            }
        }
    }
    return;
}
static void ZeroGradient(const Real Uh[restrict], Real U[restrict])
{
    for (int dim = 0; dim < DIMU; ++dim) {
        U[dim] = Uh[dim];
    }
    return;
}
/* a good practice: end file with a newline */

