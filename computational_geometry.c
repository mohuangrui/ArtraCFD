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
#include "computational_geometry.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <float.h> /* size of floating point values */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ComputeParametersAnalyticalSphere(const Space *, Polyhedron *);
static int ComputeParametersTriangulatedPolyhedron(const Space *, Polyhedron *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ComputeGeometryParameters(const Space *space, Geometry *geo)
{
    for (int m = 0; m < geo->totalM; ++m) {
        if (0 == geo->list[m].facetN) {
            ComputeParametersAnalyticalSphere(space, geo->list + m);
            continue;
        }
        ComputeParametersTriangulatedPolyhedron(space, geo->list + m);
    }
    return 0;
}
static int ComputeParametersAnalyticalSphere(const Space *space, Polyhedron *poly)
{
    const Real pi = acos(-1);
    poly->xMin = poly->xc - poly->r;
    poly->yMin = poly->yc - poly->r;
    poly->zMin = poly->zc - poly->r;
    poly->xMax = poly->xc + poly->r;
    poly->yMax = poly->yc + poly->r;
    poly->zMax = poly->zc + poly->r;
    if (COLLAPSEN != space->collapsed) { /* space dimension collapsed */
        poly->area = 2.0 * pi * poly->r; /* side area of a unit thickness cylinder */
        poly->volume = pi * poly->r * poly->r; /* volume of a unit thickness cylinder */
    } else {
        poly->area = 4.0 * pi * poly->r * poly->r; /* area of a sphere */
        poly->volume = 4.0 * pi * poly->r * poly->r * poly->r / 3.0; /* volume of a sphere */
    }
    return 0;
}
static int ComputeParametersTriangulatedPolyhedron(const Space *space, Polyhedron *poly)
{
    /* initialize parameters */
    poly->area = 0.0;
    poly->volume = 0.0;
    poly->xMin = FLT_MAX;
    poly->yMin = FLT_MAX;
    poly->zMin = FLT_MAX;
    poly->xMax = FLT_MIN;
    poly->yMax = FLT_MIN;
    poly->zMax = FLT_MIN;
    for (int n = 0; n < poly->facetN; ++n) {
        /* bounding box */
        poly->xMin = MinReal(poly->xMin, MinReal(poly->facet[n].x3, MinReal(poly->facet[n].x2, poly->facet[n].x1)));
        poly->yMin = MinReal(poly->yMin, MinReal(poly->facet[n].y3, MinReal(poly->facet[n].y2, poly->facet[n].y1)));
        poly->zMin = MinReal(poly->zMin, MinReal(poly->facet[n].z3, MinReal(poly->facet[n].z2, poly->facet[n].z1)));
        poly->xMax = MaxReal(poly->xMax, MaxReal(poly->facet[n].x3, MaxReal(poly->facet[n].x2, poly->facet[n].x1)));
        poly->yMax = MaxReal(poly->yMax, MaxReal(poly->facet[n].y3, MaxReal(poly->facet[n].y2, poly->facet[n].y1)));
        poly->zMax = MaxReal(poly->zMax, MaxReal(poly->facet[n].z3, MaxReal(poly->facet[n].z2, poly->facet[n].z1)));
        /* normal vector */
        poly->facet[n].nx = (poly->facet[n].y2 - poly->facet[n].y1) * (poly->facet[n].z3 - poly->facet[n].z1) -
            (poly->facet[n].y3 - poly->facet[n].y1) * (poly->facet[n].z2 - poly->facet[n].z1);
        poly->facet[n].ny = (poly->facet[n].x3 - poly->facet[n].x1) * (poly->facet[n].z2 - poly->facet[n].z1) -
            (poly->facet[n].x2 - poly->facet[n].x1) * (poly->facet[n].z3 - poly->facet[n].z1);
        poly->facet[n].nz = (poly->facet[n].x2 - poly->facet[n].x1) * (poly->facet[n].y3 - poly->facet[n].y1) -
            (poly->facet[n].x3 - poly->facet[n].x1) * (poly->facet[n].y2 - poly->facet[n].y1);
        poly->facet[n].s = 0.5 * sqrt(poly->facet[n].nx * poly->facet[n].nx + poly->facet[n].ny * poly->facet[n].ny + poly->facet[n].nz * poly->facet[n].nz);
        /* accumulate area and volume */
        poly->area = poly->area + poly->facet[n].s;
        poly->volume = poly->volume + poly->facet[n].x1 * poly->facet[n].nx + poly->facet[n].y1 * poly->facet[n].ny + poly->facet[n].z1 * poly->facet[n].nz;
        /* unit normal */
        poly->facet[n].nx = poly->facet[n].nx / (2.0 * poly->facet[n].s);
        poly->facet[n].ny = poly->facet[n].ny / (2.0 * poly->facet[n].s);
        poly->facet[n].nz = poly->facet[n].nz / (2.0 * poly->facet[n].s);
    }
    poly->volume = poly->volume / 6.0; /* final volume of the polyhedron */
    if (COLLAPSEN != space->collapsed) { /* space dimension collapsed */
        poly->area = poly->area - 2.0 * poly->volume; /* change to side area of a unit thickness polygon */
    }
    return 0;
}
/* a good practice: end file with a newline */

