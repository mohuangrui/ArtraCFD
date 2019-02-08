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
#include "source_term.h"
#include <string.h> /* manipulating strings */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void ComputePhi(const int tn, const int k, const int j, const int i,
        const int partn[restrict], const Node *const node,
        const Model *model, Real Phi[restrict])
{
    if (0 == model->sState) {
        memset(Phi, 0, DIMU * sizeof(*Phi));
        return;
    }
    const int idx = IndexNode(k, j, i, partn[Y], partn[X]);
    const Real *restrict U = node[idx].U[tn];
    const RealVec V = {U[1] / U[0], U[2] / U[0], U[3] / U[0]};
    const RealVec fb = {U[0] * model->g[X], U[0] * model->g[Y], U[0] * model->g[Z]};
    Phi[0] = 0.0;
    Phi[1] = fb[X];
    Phi[2] = fb[Y];
    Phi[3] = fb[Z];
    Phi[4] = Dot(fb, V);
    return;
}
/* a good practice: end file with a newline */

