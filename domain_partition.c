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
int DomainPartition(const Space *space, Partition *part)
{
    ShowInformation("Domain partitioning...");
    /*
     * Outward facing surface unit normal vector values of domain boundary, the
     * introducing of surface normal vector can provide great advantange: every
     * surface can be handled uniformly if manipulations and calculations are
     * incorporated with surface normal vector. For example, (normalX, normalY,
     * normalZ) is the outward facing unit normal vector at surface point (i,
     * j, k), then its neighbour node (ih, jh, kh) is more inner than current
     * node if (ih-i)*normalX + (jh-j)*normalY + (kh-k)*normalZ < 0.
     */
    part->normal[BCWEST][Z] = 0;
    part->normal[BCWEST][Y] = 0;
    part->normal[BCWEST][X] = -1;

    part->normal[BCEAST][Z] = 0;
    part->normal[BCEAST][Y] = 0;
    part->normal[BCEAST][X] = 1;

    part->normal[BCSOUTH][Z] = 0;
    part->normal[BCSOUTH][Y] = -1;
    part->normal[BCSOUTH][X] = 0;

    part->normal[BCNORTH][Z] = 0;
    part->normal[BCNORTH][Y] = 1;
    part->normal[BCNORTH][X] = 0;

    part->normal[BCFRONT][Z] = -1;
    part->normal[BCFRONT][Y] = 0;
    part->normal[BCFRONT][X] = 0;

    part->normal[BCBACK][Z] = 1;
    part->normal[BCBACK][Y] = 0;
    part->normal[BCBACK][X] = 0;
    /*
     * Assign values to each index control array of each partition.
     *
     * (In this program, min and max are used for domain index identifying,
     * therefore they are widely used in for loop control. To reduce the for loop
     * operations, we always use a reachable value for min, and unreachable value
     * for the max. To count the total valid objects, simply do max - min).
     *
     * Entire domain (includes exterior ghost): min = 0; max = m + 2 * ng;
     * Lower exterior ghost: min = 0; max = ng;
     * Normal nodes: min = ng; max = m + ng;
     *      Lower boundary: min = ng; max = ng + 1;
     *      Interior cells(node layers): min = ng + 1; max = m + ng - 1;
     *      Upper Boundary: min = m + ng - 1; max = m + ng;
     * Upper exterior ghost: min = m + ng; max = n;
     *
     * Note that for each direction, its boundary nodes and exterior ghost
     * nodes will only extent out from the interior cells at that direction
     * and will not extent on other directions, that is, they form cross
     * like shapes in space without corner parts.
     *
     * Apparently it represents a certain problem, since it is not quite clear
     * how to set corner values (if there is no adjacent grid block). The values
     * are not required by the standard cross-type discretisation stencil.
     * However, they may become necessary for the computation of gradients
     * (viscous fluxes), or for transfer operators within multigrid.
     * Usually, an averaging of the values from the adjacent “regular” dummy
     * cells is sufficient. However, it's really complicated for 3D domian with
     * multiple ghost layers: there are 12 edge corner blocks and 8 vertex
     * corner block need to be handled. This problem is left here.
     *
     */
    part->n[PIN][Z][MIN] = space->ng + 1;
    part->n[PIN][Z][MAX] = space->m[Z] + space->ng - 1;
    part->n[PIN][Y][MIN] = space->ng + 1;
    part->n[PIN][Y][MAX] = space->m[Y] + space->ng - 1;
    part->n[PIN][X][MIN] = space->ng + 1;
    part->n[PIN][X][MAX] = space->m[X] + space->ng - 1;

    part->n[PWB][Z][MIN] = space->ng + 1;
    part->n[PWB][Z][MAX] = space->m[Z] + space->ng - 1;
    part->n[PWB][Y][MIN] = space->ng + 1;
    part->n[PWB][Y][MAX] = space->m[Y] + space->ng - 1;
    part->n[PWB][X][MIN] = space->ng;
    part->n[PWB][X][MAX] = space->ng + 1;

    part->n[PEB][Z][MIN] = space->ng + 1;
    part->n[PEB][Z][MAX] = space->m[Z] + space->ng - 1;
    part->n[PEB][Y][MIN] = space->ng + 1;
    part->n[PEB][Y][MAX] = space->m[Y] + space->ng - 1;
    part->n[PEB][X][MIN] = space->m[X] + space->ng - 1;
    part->n[PEB][X][MAX] = space->m[X] + space->ng;

    part->n[PSB][Z][MIN] = space->ng + 1;
    part->n[PSB][Z][MAX] = space->m[Z] + space->ng - 1;
    part->n[PSB][Y][MIN] = space->ng;
    part->n[PSB][Y][MAX] = space->ng + 1;
    part->n[PSB][X][MIN] = space->ng + 1;
    part->n[PSB][X][MAX] = space->m[X] + space->ng - 1;

    part->n[PNB][Z][MIN] = space->ng + 1;
    part->n[PNB][Z][MAX] = space->m[Z] + space->ng - 1;
    part->n[PNB][Y][MIN] = space->m[Y] + space->ng - 1;
    part->n[PNB][Y][MAX] = space->m[Y] + space->ng;
    part->n[PNB][X][MIN] = space->ng + 1;
    part->n[PNB][X][MAX] = space->m[X] + space->ng - 1;

    part->n[PFB][Z][MIN] = space->ng;
    part->n[PFB][Z][MAX] = space->ng + 1;
    part->n[PFB][Y][MIN] = space->ng + 1;
    part->n[PFB][Y][MAX] = space->m[Y] + space->ng - 1;
    part->n[PFB][X][MIN] = space->ng + 1;
    part->n[PFB][X][MAX] = space->m[X] + space->ng - 1;

    part->n[PBB][Z][MIN] = space->m[Z] + space->ng - 1;
    part->n[PBB][Z][MAX] = space->m[Z] + space->ng;
    part->n[PBB][Y][MIN] = space->ng + 1;
    part->n[PBB][Y][MAX] = space->m[Y] + space->ng - 1;
    part->n[PBB][X][MIN] = space->ng + 1;
    part->n[PBB][X][MAX] = space->m[X] + space->ng - 1;

    part->n[PWG][Z][MIN] = space->ng + 1;
    part->n[PWG][Z][MAX] = space->m[Z] + space->ng - 1;
    part->n[PWG][Y][MIN] = space->ng + 1;
    part->n[PWG][Y][MAX] = space->m[Y] + space->ng - 1;
    part->n[PWG][X][MIN] = 0;
    part->n[PWG][X][MAX] = space->ng;

    part->n[PEG][Z][MIN] = space->ng + 1;
    part->n[PEG][Z][MAX] = space->m[Z] + space->ng - 1;
    part->n[PEG][Y][MIN] = space->ng + 1;
    part->n[PEG][Y][MAX] = space->m[Y] + space->ng - 1;
    part->n[PEG][X][MIN] = space->m[X] + space->ng;
    part->n[PEG][X][MAX] = space->m[X] + 2 * space->ng;

    part->n[PSG][Z][MIN] = space->ng + 1;
    part->n[PSG][Z][MAX] = space->m[Z] + space->ng - 1;
    part->n[PSG][Y][MIN] = 0;
    part->n[PSG][Y][MAX] = space->ng;
    part->n[PSG][X][MIN] = space->ng + 1;
    part->n[PSG][X][MAX] = space->m[X] + space->ng - 1;

    part->n[PNG][Z][MIN] = space->ng + 1;
    part->n[PNG][Z][MAX] = space->m[Z] + space->ng - 1;
    part->n[PNG][Y][MIN] = space->m[Y] + space->ng;
    part->n[PNG][Y][MAX] = space->m[Y] + 2 * space->ng;
    part->n[PNG][X][MIN] = space->ng + 1;
    part->n[PNG][X][MAX] = space->m[X] + space->ng - 1;

    part->n[PFG][Z][MIN] = 0;
    part->n[PFG][Z][MAX] = space->ng;
    part->n[PFG][Y][MIN] = space->ng + 1;
    part->n[PFG][Y][MAX] = space->m[Y] + space->ng - 1;
    part->n[PFG][X][MIN] = space->ng + 1;
    part->n[PFG][X][MAX] = space->m[X] + space->ng - 1;

    part->n[PBG][Z][MIN] = space->m[Z] + space->ng;
    part->n[PBG][Z][MAX] = space->m[Z] + 2 * space->ng;
    part->n[PBG][Y][MIN] = space->ng + 1;
    part->n[PBG][Y][MAX] = space->m[Y] + space->ng - 1;
    part->n[PBG][X][MIN] = space->ng + 1;
    part->n[PBG][X][MAX] = space->m[X] + space->ng - 1;
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

