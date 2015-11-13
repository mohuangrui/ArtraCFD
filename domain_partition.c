/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) 2014-2018 Huangrui Mo <huangrui.mo@gmail.com>              *
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
     * Assign values to each index control array of each inner partition.
     *
     * Member nx, ny, nz in Space is the total number of node layers.
     * The detailed node range of the domain will be the following:
     *
     * (In this program, Sub and Sup are used for domain range identifying,
     * therefore they are widely used in for loop control. To reduce the for loop
     * operations, we always use a reachable value for Sub, and unreachable value
     * for the Sup. To count the total valid objects, simply do Sup - Sub).
     *
     * Entire domain (includes exterior ghost): Sub = 0; Sup = n + 2 * ng;
     * Lower exterior ghost: Sub = 0; Sup = ng;
     * Normal nodes: Sub = ng; Sup = n + ng;
     *      Lower boundary: Sub = ng; Sup = ng + 1;
     *      Interior cells(node layers): Sub = ng + 1; Sup = n + ng - 1;
     *      Upper Boundary: Sub = n + ng - 1; Sup = n + ng;
     * Upper exterior ghost: Sub = n + ng; Sup = n + 2 *ng;
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
    part->range[PIN][KSUB] = space->ng + 1;
    part->range[PIN][KSUP] = space->nz + space->ng - 1;
    part->range[PIN][JSUB] = space->ng + 1;
    part->range[PIN][JSUP] = space->ny + space->ng - 1;
    part->range[PIN][ISUB] = space->ng + 1;
    part->range[PIN][ISUP] = space->nx + space->ng - 1;

    part->range[PWB][KSUB] = space->ng + 1;
    part->range[PWB][KSUP] = space->nz + space->ng - 1;
    part->range[PWB][JSUB] = space->ng + 1;
    part->range[PWB][JSUP] = space->ny + space->ng - 1;
    part->range[PWB][ISUB] = space->ng;
    part->range[PWB][ISUP] = space->ng + 1;

    part->range[PEB][KSUB] = space->ng + 1;
    part->range[PEB][KSUP] = space->nz + space->ng - 1;
    part->range[PEB][JSUB] = space->ng + 1;
    part->range[PEB][JSUP] = space->ny + space->ng - 1;
    part->range[PEB][ISUB] = space->nx + space->ng - 1;
    part->range[PEB][ISUP] = space->nx + space->ng;

    part->range[PSB][KSUB] = space->ng + 1;
    part->range[PSB][KSUP] = space->nz + space->ng - 1;
    part->range[PSB][JSUB] = space->ng;
    part->range[PSB][JSUP] = space->ng + 1;
    part->range[PSB][ISUB] = space->ng + 1;
    part->range[PSB][ISUP] = space->nx + space->ng - 1;

    part->range[PNB][KSUB] = space->ng + 1;
    part->range[PNB][KSUP] = space->nz + space->ng - 1;
    part->range[PNB][JSUB] = space->ny + space->ng - 1;
    part->range[PNB][JSUP] = space->ny + space->ng;
    part->range[PNB][ISUB] = space->ng + 1;
    part->range[PNB][ISUP] = space->nx + space->ng - 1;

    part->range[PFB][KSUB] = space->ng;
    part->range[PFB][KSUP] = space->ng + 1;
    part->range[PFB][JSUB] = space->ng + 1;
    part->range[PFB][JSUP] = space->ny + space->ng - 1;
    part->range[PFB][ISUB] = space->ng + 1;
    part->range[PFB][ISUP] = space->nx + space->ng - 1;

    part->range[PBB][KSUB] = space->nz + space->ng - 1;
    part->range[PBB][KSUP] = space->nz + space->ng;
    part->range[PBB][JSUB] = space->ng + 1;
    part->range[PBB][JSUP] = space->ny + space->ng - 1;
    part->range[PBB][ISUB] = space->ng + 1;
    part->range[PBB][ISUP] = space->nx + space->ng - 1;

    part->range[PWG][KSUB] = space->ng + 1;
    part->range[PWG][KSUP] = space->nz + space->ng - 1;
    part->range[PWG][JSUB] = space->ng + 1;
    part->range[PWG][JSUP] = space->ny + space->ng - 1;
    part->range[PWG][ISUB] = 0;
    part->range[PWG][ISUP] = space->ng;

    part->range[PEG][KSUB] = space->ng + 1;
    part->range[PEG][KSUP] = space->nz + space->ng - 1;
    part->range[PEG][JSUB] = space->ng + 1;
    part->range[PEG][JSUP] = space->ny + space->ng - 1;
    part->range[PEG][ISUB] = space->nx + space->ng;
    part->range[PEG][ISUP] = space->nx + 2 * space->ng;

    part->range[PSG][KSUB] = space->ng + 1;
    part->range[PSG][KSUP] = space->nz + space->ng - 1;
    part->range[PSG][JSUB] = 0;
    part->range[PSG][JSUP] = space->ng;
    part->range[PSG][ISUB] = space->ng + 1;
    part->range[PSG][ISUP] = space->nx + space->ng - 1;

    part->range[PNG][KSUB] = space->ng + 1;
    part->range[PNG][KSUP] = space->nz + space->ng - 1;
    part->range[PNG][JSUB] = space->ny + space->ng;
    part->range[PNG][JSUP] = space->ny + 2 * space->ng;
    part->range[PNG][ISUB] = space->ng + 1;
    part->range[PNG][ISUP] = space->nx + space->ng - 1;

    part->range[PFG][KSUB] = 0;
    part->range[PFG][KSUP] = space->ng;
    part->range[PFG][JSUB] = space->ng + 1;
    part->range[PFG][JSUP] = space->ny + space->ng - 1;
    part->range[PFG][ISUB] = space->ng + 1;
    part->range[PFG][ISUP] = space->nx + space->ng - 1;

    part->range[PBG][KSUB] = space->nz + space->ng;
    part->range[PBG][KSUP] = space->nz + 2 * space->ng;
    part->range[PBG][JSUB] = space->ng + 1;
    part->range[PBG][JSUP] = space->ny + space->ng - 1;
    part->range[PBG][ISUB] = space->ng + 1;
    part->range[PBG][ISUP] = space->nx + space->ng - 1;
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

