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
int DomainPartition(Space *space)
{
    Partition *part = &(space->part);
    /*
     * Outward facing surface unit normal vector values of domain boundary, the
     * introducing of surface normal vector can provide great advantange: every
     * surface can be handled uniformly if manipulations and calculations are
     * incorporated with surface normal vector. For example, (NX, NY, NZ) is the
     * outward facing unit normal vector at surface point (i, j, k), then its
     * neighbour node (ih, jh, kh) is more inner than current node if 
     * (ih - i) * NX + (jh - j) * NY + (kh - k) * NZ < 0.
     */
    part->N[PWB][X] = -1;
    part->N[PWB][Y] = 0;
    part->N[PWB][Z] = 0;

    part->N[PEB][X] = 1;
    part->N[PEB][Y] = 0;
    part->N[PEB][Z] = 0;

    part->N[PSB][X] = 0;
    part->N[PSB][Y] = -1;
    part->N[PSB][Z] = 0;

    part->N[PNB][X] = 0;
    part->N[PNB][Y] = 1;
    part->N[PNB][Z] = 0;

    part->N[PFB][X] = 0;
    part->N[PFB][Y] = 0;
    part->N[PFB][Z] = -1;

    part->N[PBB][X] = 0;
    part->N[PBB][Y] = 0;
    part->N[PBB][Z] = 1;
    /*
     * Assign values to each index control array of each partition.
     *
     * (In this program, min and max are used for domain index identifying,
     * therefore they are widely used in for loop control. To reduce the for loop
     * operations, we always use a reachable value for min, and unreachable value
     * for the max. To count the total valid objects, simply do max - min).
     *
     * Index range
     *
     * Entire domain (includes exterior ghost): [0, n);
     * Lower exterior ghost: [0, ng);
     * Normal nodes: [ng, n-ng);
     *      Lower boundary: [ng, ng+1);
     *      Interior node layers: [ng+1, n-ng-1);
     *      Upper Boundary: [n-ng-1, n-ng);
     * Upper exterior ghost: [n-ng, n);
     *
     * For a 3D domain with multiple ghost layers, there are 12 edge corner
     * blocks and 8 vertex corner blocks formed by ghost nodes. It is not
     * clear how to set corner ghost values (if no adjacent grid block).
     * The values are not required by cross-type discretization stencils.
     * However, they may become necessary for the computation of mixed 
     * derivatives as in viscous fluxes, or for transfer operators within
     * multigrid. Usually, an averaging of the values from the adjacent
     * cells is sufficient. However, it's really complicated in 3D space.
     *
     * When mixed derivatives are discretized by second order schemes, 
     * only the nearest corner nodes are required for each computational 
     * node. Therefore, the global exterior ghost nodes will only extent
     * out from the interior nodes at that direction and will not extent
     * on other parts, that is, they form cross like shapes in space 
     * without corner parts. For the global boundary nodes, they form a 
     * square with corner nodes to provide stencils for mix derivatives.
     */
    part->ns[PIN][X][MIN] = part->ng + 1;
    part->ns[PIN][X][MAX] = part->n[X] - part->ng - 1;
    part->ns[PIN][Y][MIN] = part->ng + 1;
    part->ns[PIN][Y][MAX] = part->n[Y] - part->ng - 1;
    part->ns[PIN][Z][MIN] = part->ng + 1;
    part->ns[PIN][Z][MAX] = part->n[Z] - part->ng - 1;

    part->ns[PWB][X][MIN] = part->ng;
    part->ns[PWB][X][MAX] = part->ng + 1;
    part->ns[PWB][Y][MIN] = part->ng;
    part->ns[PWB][Y][MAX] = part->n[Y] - part->ng;
    part->ns[PWB][Z][MIN] = part->ng;
    part->ns[PWB][Z][MAX] = part->n[Z] - part->ng;

    part->ns[PEB][X][MIN] = part->n[X] - part->ng - 1;
    part->ns[PEB][X][MAX] = part->n[X] - part->ng;
    part->ns[PEB][Y][MIN] = part->ng;
    part->ns[PEB][Y][MAX] = part->n[Y] - part->ng;
    part->ns[PEB][Z][MIN] = part->ng;
    part->ns[PEB][Z][MAX] = part->n[Z] - part->ng;

    part->ns[PSB][X][MIN] = part->ng;
    part->ns[PSB][X][MAX] = part->n[X] - part->ng;
    part->ns[PSB][Y][MIN] = part->ng;
    part->ns[PSB][Y][MAX] = part->ng + 1;
    part->ns[PSB][Z][MIN] = part->ng;
    part->ns[PSB][Z][MAX] = part->n[Z] - part->ng;

    part->ns[PNB][X][MIN] = part->ng;
    part->ns[PNB][X][MAX] = part->n[X] - part->ng;
    part->ns[PNB][Y][MIN] = part->n[Y] - part->ng - 1;
    part->ns[PNB][Y][MAX] = part->n[Y] - part->ng;
    part->ns[PNB][Z][MIN] = part->ng;
    part->ns[PNB][Z][MAX] = part->n[Z] - part->ng;

    part->ns[PFB][X][MIN] = part->ng;
    part->ns[PFB][X][MAX] = part->n[X] - part->ng;
    part->ns[PFB][Y][MIN] = part->ng;
    part->ns[PFB][Y][MAX] = part->n[Y] - part->ng;
    part->ns[PFB][Z][MIN] = part->ng;
    part->ns[PFB][Z][MAX] = part->ng + 1;

    part->ns[PBB][X][MIN] = part->ng;
    part->ns[PBB][X][MAX] = part->n[X] - part->ng;
    part->ns[PBB][Y][MIN] = part->ng;
    part->ns[PBB][Y][MAX] = part->n[Y] - part->ng;
    part->ns[PBB][Z][MIN] = part->n[Z] - part->ng - 1;
    part->ns[PBB][Z][MAX] = part->n[Z] - part->ng;

    part->ns[PWG][X][MIN] = 0;
    part->ns[PWG][X][MAX] = part->ng;
    part->ns[PWG][Y][MIN] = part->ng + 1;
    part->ns[PWG][Y][MAX] = part->n[Y] - part->ng - 1;
    part->ns[PWG][Z][MIN] = part->ng + 1;
    part->ns[PWG][Z][MAX] = part->n[Z] - part->ng - 1;

    part->ns[PEG][X][MIN] = part->n[X] - part->ng;
    part->ns[PEG][X][MAX] = part->n[X];
    part->ns[PEG][Y][MIN] = part->ng + 1;
    part->ns[PEG][Y][MAX] = part->n[Y] - part->ng - 1;
    part->ns[PEG][Z][MIN] = part->ng + 1;
    part->ns[PEG][Z][MAX] = part->n[Z] - part->ng - 1;

    part->ns[PSG][X][MIN] = part->ng + 1;
    part->ns[PSG][X][MAX] = part->n[X] - part->ng - 1;
    part->ns[PSG][Y][MIN] = 0;
    part->ns[PSG][Y][MAX] = part->ng;
    part->ns[PSG][Z][MIN] = part->ng + 1;
    part->ns[PSG][Z][MAX] = part->n[Z] - part->ng - 1;

    part->ns[PNG][X][MIN] = part->ng + 1;
    part->ns[PNG][X][MAX] = part->n[X] - part->ng - 1;
    part->ns[PNG][Y][MIN] = part->n[Y] - part->ng;
    part->ns[PNG][Y][MAX] = part->n[Y];
    part->ns[PNG][Z][MIN] = part->ng + 1;
    part->ns[PNG][Z][MAX] = part->n[Z] - part->ng - 1;

    part->ns[PFG][X][MIN] = part->ng + 1;
    part->ns[PFG][X][MAX] = part->n[X] - part->ng - 1;
    part->ns[PFG][Y][MIN] = part->ng + 1;
    part->ns[PFG][Y][MAX] = part->n[Y] - part->ng - 1;
    part->ns[PFG][Z][MIN] = 0;
    part->ns[PFG][Z][MAX] = part->ng;

    part->ns[PBG][X][MIN] = part->ng + 1;
    part->ns[PBG][X][MAX] = part->n[X] - part->ng - 1;
    part->ns[PBG][Y][MIN] = part->ng + 1;
    part->ns[PBG][Y][MAX] = part->n[Y] - part->ng - 1;
    part->ns[PBG][Z][MIN] = part->n[Z] - part->ng;
    part->ns[PBG][Z][MAX] = part->n[Z];

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
    return 0;
}
/* a good practice: end file with a newline */

