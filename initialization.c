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
#include "initialization.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "boundary_treatment.h"
#include "data_stream.h"
#include "geometry_stream.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NonRestartInitializer(Real *U, const Space *, const Geometry *,
        const Partition *, const Flow *);
static int ApplyRegionalInitializer(const int, Real *U, const Space *, 
        const Partition *, const Flow *);
static int RestartInitializer(Real *U, const Space *, const Geometry *, Time *, 
        const Partition *, const Flow *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function initializes the entire flow field. Initialization will be 
 * done differently determined by the restart status.
 */
int InitializeFlowField(Real *U, const Space *space, const Geometry *geometry,
        Time *time, const Partition *part, const Flow *flow)
{
    ShowInformation("Initializing flow field...");
    if (0 == time->restart) { /* non restart */
        NonRestartInitializer(U, space, geometry, part, flow);
        /* if this is a first run, output initial data */
        WriteComputedData(U, space, time, part, flow);
        WriteGeometryData(geometry, time);
    } else {
        RestartInitializer(U, space, geometry, time, part, flow);
    }
    ShowInformation("Session End");
    return 0;
}
/*
 * The non-restart initialization will assign values to field variables.
 */
static int NonRestartInitializer(Real *U, const Space *space, const Geometry *geometry, 
        const Partition *part, const Flow *flow)
{
    ShowInformation("  Non-restart run initializing...");
    /* extract initial values */
    const Real Uo[DIMUo] = {
        part->valueBC[0][0],
        part->valueBC[0][1],
        part->valueBC[0][2],
        part->valueBC[0][3],
        part->valueBC[0][4]};
    /*
     * Initialize the interior field
     */
    int idx = 0; /* linear array index math variable */
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                idx = IndexMath(k, j, i, space) * DIMU;
                ConservativeByPrimitive(U, idx, Uo, flow);
            }
        }
    }
    /*
     * Regional initializer for specific flow regions
     */
    for (int n = 0; n < part->tallyIC; ++n) {
        ApplyRegionalInitializer(n, U, space, part, flow);
    }
    /*
     * Boundary conditions and treatments to obtain an entire initialized flow field
     */
    BoundaryCondtionsAndTreatments(U, space, geometry, part, flow);
    return 0;
}
/*
 * The handling of regional initialization for specific flow region is achieved
 * through the cooperation of three data structures:
 * The tallyIC counts the number of initializers.
 * The typeIC array keeps a list of the types of regional initialization.
 * The valueIC array stored the information of the corresponding IC type.
 * IC types and corresponding information entries:
 * 1: plane (x, y, z, normalX, normalY, normalZ, rho, u, v, w, p)
 * 2: sphere (x, y, z, r, rho, u, v, w, p)
 * 3: box (xmin, ymin, zmin, xmax, ymax, zmax, rho, u, v, w, p)
 * 4: cylinder (x1, y1, z1, x2, y2, z2, r, rho, u, v, w, p)
 */
static int ApplyRegionalInitializer(const int n, Real *U, const Space *space, 
        const Partition *part, const Flow *flow)
{
    /*
     * Acquire the specialized information data entries
     */
    /* the fix index part */
    const Real x1 = part->valueIC[n][0];
    const Real y1 = part->valueIC[n][1];
    const Real z1 = part->valueIC[n][2];
    const Real Uo[DIMUo] = {
        part->valueIC[n][ENTRYIC-5],
        part->valueIC[n][ENTRYIC-4],
        part->valueIC[n][ENTRYIC-3],
        part->valueIC[n][ENTRYIC-2],
        part->valueIC[n][ENTRYIC-1]};
    /* the vary part */
    Real r = 0.0;
    Real x2 = 0.0;
    Real y2 = 0.0;
    Real z2 = 0.0;
    Real normalZ = 0.0;
    Real normalY = 0.0;
    Real normalX = 0.0;
    switch (part->typeIC[n]) {
        case 1: /* plane */
            normalX = part->valueIC[n][3];
            normalY = part->valueIC[n][4];
            normalZ = part->valueIC[n][5];
            break;
        case 2: /* sphere */
            r = part->valueIC[n][3];
            break;
        case 3: /* box */
            x2 = part->valueIC[n][3];
            y2 = part->valueIC[n][4];
            z2 = part->valueIC[n][5];
            break;
        case 4: /* cylinder */
            x2 = part->valueIC[n][3];
            y2 = part->valueIC[n][4];
            z2 = part->valueIC[n][5];
            r = part->valueIC[n][6];
            break;
        default:
            break;
    }
    /*
     * Apply initial values for nodes that meets condition
     */
    int idx = 0;
    int flag = 0; /* control flag for whether current node in the region */
    const Real xVec = x2 - x1;
    const Real yVec = y2 - y1;
    const Real zVec = z2 - z1;
    const Real lVecSquare = xVec * xVec + yVec * yVec + zVec * zVec;
    Real x = 0.0;
    Real y = 0.0;
    Real z = 0.0;
    Real xh = 0.0;
    Real yh = 0.0;
    Real zh = 0.0;
    Real dot = 0.0;
    Real distSquare = 0.0;
    for (int k = part->kSub[0]; k < part->kSup[0]; ++k) {
        for (int j = part->jSub[0]; j < part->jSup[0]; ++j) {
            for (int i = part->iSub[0]; i < part->iSup[0]; ++i) {
                flag = 0; /* always initialize flag to zero */
                x = ComputeX(i, space);
                y = ComputeY(j, space);
                z = ComputeZ(k, space);
                xh = (x - x1);
                yh = (y - y1);
                zh = (z - z1);
                switch (part->typeIC[n]) {
                    case 1: /* plane */
                        dot = xh * normalX + yh * normalY + zh * normalZ;
                        if (0 <= dot) { /* on the normal direction or the plane */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case 2: /* sphere */
                        if (0 >= (xh * xh + yh * yh + zh * zh - r * r)) { /* in or on the sphere */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case 3: /* box */
                        xh = xh * (x - x2);
                        yh = yh * (y - y2);
                        zh = zh * (z - z2);
                        if ((0 >= xh) && (0 >= yh) && (0 >= zh)) { /* in or on the box */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    case 4: /* cylinder */
                        dot = xh * xVec + yh * yVec + zh * zVec;
                        if ((0 > dot) || (lVecSquare < dot)) { /* outside the two ends */
                            break;
                        }
                        distSquare = xh * xh + yh * yh + zh * zh - dot * dot / lVecSquare;
                        if (r * r >= distSquare) { /* in or on the cylinder */
                            flag = 1; /* set flag to true */
                        }
                        break;
                    default:
                        break;
                }
                if (1 == flag) { /* current node meets the condition */
                    idx = IndexMath(k, j, i, space) * DIMU;
                    ConservativeByPrimitive(U, idx, Uo, flow);
                }
            }
        }
    }
    return 0;
}
/*
 * If this is a restart run, then initialize flow field by reading field data
 * from restart files.
 */
static int RestartInitializer(Real *U, const Space *space, const Geometry *geometry, 
        Time *time, const Partition *part, const Flow *flow)
{
    ShowInformation("  Restart run initializing...");
    /*
     * Load data from Ensight restart files.
     */
    LoadComputedData(U, space, time, part, flow);
    /*
     * Boundary conditions and treatments to obtain an entire initialized flow field
     */
    BoundaryCondtionsAndTreatments(U, space, geometry, part, flow);
    return 0;
}
/* a good practice: end file with a newline */

