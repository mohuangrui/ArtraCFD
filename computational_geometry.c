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
#include "computational_geometry.h"
#include <stdio.h> /* standard library for input and output */
#include <math.h> /* common mathematical functions */
#include <float.h> /* size of floating point values */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ComputeParametersAnalyticalSphere(Polyhedron *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ComputeGeometryParameters(Geometry *geo)
{
    for (int m = 0; m < geo->totalM; ++m) {
        if (0 == geo->list[m].facetN) {
            ComputeParametersAnalyticalSphere(geo->list + m);
            continue;
        }
    }
}
static int ComputeParametersAnalyticalSphere(Polyhedron *poly)
{
    const Real pi = acos(-1);
    poly->xMin = poly->xc - poly->r;
    poly->yMin = poly->yc - poly->r;
    poly->zMin = poly->zc - poly->r;
    poly->xMax = poly->xc + poly->r;
    poly->yMax = poly->yc + poly->r;
    poly->zMax = poly->zc + poly->r;
    poly->area = 4.0 * pi * poly->r * poly->r;
    poly->volume = poly->area * poly->r / 3.0;
    return 0;
}
static int ComputeParametersTriangulatedPolyhedron(Polyhedron *poly)
{
    Real magNormal = 0.0;
    /* initialize parameters */
    poly->area = 0.0;
    poly->volume = 0.0;
    poly->xMin = FLT_MIN;
    poly->yMin = FLT_MIN;
    poly->zMin = FLT_MIN;
    poly->xMax = FLT_MAX;
    poly->yMax = FLT_MAX;
    poly->zMax = FLT_MAX;
    for (int n = 0; n < poly->facetN; ++n) {
        /* normal vector */
        poly->facet[n].nx = (poly->facet[n].y2 - poly->facet[n].y1) * (poly->facet[n].z3 - poly->facet[n].z1) -
            (poly->facet[n].y3 - poly->facet[n].y1) * (poly->facet[n].z2 - poly->facet[n].z1);
        poly->facet[n].ny = (poly->facet[n].x3 - poly->facet[n].x1) * (poly->facet[n].z2 - poly->facet[n].z1) -
            (poly->facet[n].x2 - poly->facet[n].x1) * (poly->facet[n].z3 - poly->facet[n].z1);
        poly->facet[n].nz = (poly->facet[n].x2 - poly->facet[n].x1) * (poly->facet[n].y3 - poly->facet[n].y1) -
            (poly->facet[n].x3 - poly->facet[n].x1) * (poly->facet[n].y2 - poly->facet[n].y1);
        magNormal = sqrt(poly->facet[n].nx * poly->facet[n].nx + poly->facet[n].ny * poly->facet[n].ny + poly->facet[n].nz * poly->facet[n].nz);
        /* accumulate area */
        poly->area = poly->area + 0.5 * magNormal;
    }
    return 0;
}
/* a good practice: end file with a newline */

