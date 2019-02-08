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
#include "domain_partition.h"
#include "commons.h"
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void PartitionDomain(Space *space)
{
    Partition *const part = &(space->part);
    /*
     * Outward facing surface unit normal vector of domain boundary
     * Surface normal vector can provide great advantage: every surface can
     * be treated uniformly if operations are incorporated with surface
     * normal vector. For example, (NX, NY, NZ) is the outward facing unit
     * normal vector at surface point (i, j, k), then its neighbour node
     * (ih, jh, kh) is more inner than current node if
     * (ih - i) * NX + (jh - j) * NY + (kh - k) * NZ < 0.
     */
    const int N[NBC][DIMS] = {{0, 0, 0}, {-1, 0, 0}, {1, 0, 0},
        {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};
    for (int p = 0; p < NBC; ++p) {
        for (int s = 0; s < DIMS; ++s) {
            part->N[p][s] = N[p][s];
        }
    }
    /*
     * Define the index control array of each partition
     * min and max are used to define domain index range and are widely
     * used in for loop control. To reduce loop comparisons, use a
     * reachable value for min, and unreachable value for max. To count
     * the total valid objects, simply do (max - min).
     *
     * Index range
     * Entire domain (includes exterior ghost): [0, n);
     * Lower exterior ghost: [0, ng);
     * Normal nodes: [ng, n-ng);
     *      Lower boundary: [ng, ng+1);
     *      Interior node layers: [ng+1, n-ng-1);
     *      Upper Boundary: [n-ng-1, n-ng);
     * Upper exterior ghost: [n-ng, n);
     *
     * The computational domain are categorized into three main regions:
     * interior region, boundary region, and ghost region. In general,
     * the interior region is where normal computation is performed;
     * the boundary region is where physical boundary condition is enforced;
     * the ghost region is where numerical boundary treatment is performed.
     *
     * When the periodic boundary condition is applied to a direction, then
     * the interior region in this direction should be extended to include
     * the boundary region, which will then participate normal computation.
     */
    for (int s = 0, q = PWB; s < DIMS; ++s, q = q + 2) {
        /* interior region */
        if (PERIODIC == part->typeBC[q]) { /* extended interior range */
            part->ns[PIN][s][MIN] = part->ng[s];
            part->ns[PIN][s][MAX] = part->n[s] - part->ng[s];
        } else { /* normal interior range */
            part->ns[PIN][s][MIN] = part->ng[s] + 1;
            part->ns[PIN][s][MAX] = part->n[s] - part->ng[s] - 1;
        }
        /* boundary box */
        for (int p = PWB; p <= PBB; ++p) {
            part->ns[p][s][MIN] = part->ng[s];
            part->ns[p][s][MAX] = part->n[s] - part->ng[s];
        }
        /* ghost box */
        for (int p = PWG; p <= PBG; ++p) {
            part->ns[p][s][MIN] = 0;
            part->ns[p][s][MAX] = part->n[s];
        }
        /* physical region */
        part->ns[PHY][s][MIN] = part->ng[s];
        part->ns[PHY][s][MAX] = part->n[s] - part->ng[s];
        /* all region */
        part->ns[PAL][s][MIN] = 0;
        part->ns[PAL][s][MAX] = part->n[s];
    }
    /* adjust boundary box to the boundary region of each direction */
    for (int p = PWB, s = 0; p <= PBB; p = p + 2, ++s) {
        part->ns[p][s][MAX] = part->ns[p][s][MIN] + 1;
    }
    for (int p = PEB, s = 0; p <= PBB; p = p + 2, ++s) {
        part->ns[p][s][MIN] = part->ns[p][s][MAX] - 1;
    }
    /* adjust ghost box to the ghost region of each direction */
    for (int p = PWG, s = 0; p <= PBG; p = p + 2, ++s) {
        part->ns[p][s][MAX] = part->ng[s];
    }
    for (int p = PEG, s = 0; p <= PBG; p = p + 2, ++s) {
        part->ns[p][s][MIN] = part->n[s] - part->ng[s];
    }
    /* computational node range with dimension priority */
    const int np[DIMS][DIMS][LIMIT] = {
        {
            {part->ns[PIN][X][MIN], part->ns[PIN][X][MAX]},
            {part->ns[PIN][Y][MIN], part->ns[PIN][Y][MAX]},
            {part->ns[PIN][Z][MIN], part->ns[PIN][Z][MAX]}
        },
        {
            {part->ns[PIN][Y][MIN], part->ns[PIN][Y][MAX]},
            {part->ns[PIN][X][MIN], part->ns[PIN][X][MAX]},
            {part->ns[PIN][Z][MIN], part->ns[PIN][Z][MAX]}
        },
        {
            {part->ns[PIN][Z][MIN], part->ns[PIN][Z][MAX]},
            {part->ns[PIN][X][MIN], part->ns[PIN][X][MAX]},
            {part->ns[PIN][Y][MIN], part->ns[PIN][Y][MAX]}
        }
    };
    for (int n = 0; n < DIMS; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            for (int m = 0; m < LIMIT; ++m) {
                part->np[n][s][m] = np[n][s][m];
            }
        }
    }
    /* search path for interfacial node */
    const int path[PATHN][DIMS] = { /* searching path */
        {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1},
        {-1, -1, 0}, {1, -1, 0}, {-1, 1, 0}, {1, 1, 0},
        {-1, 0, -1}, {1, 0, -1}, {-1, 0, 1}, {1, 0, 1},
        {0, -1, -1}, {0, 1, -1}, {0, -1, 1}, {0, 1, 1},
        {-2, 0, 0}, {2, 0, 0}, {0, -2, 0}, {0, 2, 0}, {0, 0, -2}, {0, 0, 2},
        {-3, 0, 0}, {3, 0, 0}, {0, -3, 0}, {0, 3, 0}, {0, 0, -3}, {0, 0, 3}
    };
    for (int n = 0; n < PATHN; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            part->path[n][s] = path[n][s];
        }
    }
    const int base = 6; /* base directions in searching path */
    part->pathSep[1] = base; /* end index for layer 1 */
    part->pathSep[2] = part->pathSep[1] + 18; /* end index for layer 2 */
    part->pathSep[3] = part->pathSep[2] + base; /* end index for layer 3 */
    /* max search path for a spatial scheme */
    part->pathSep[0] = part->pathSep[part->gl];
    return;
}
/* a good practice: end file with a newline */

