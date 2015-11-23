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
    for (int s = 0; s < DIMS; ++s) {
        poly->box[s][MIN] = poly->O[s] - poly->r;
        poly->box[s][MAX] = poly->O[s] + poly->r;
    }
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
    RealVector P1P2 = {0.0};
    RealVector P1P3 = {0.0};
    poly->area = 0.0;
    poly->volume = 0.0;
    for (int s = 0; s < DIMS; ++s) {
        poly->box[s][MIN] = FLT_MAX;
        poly->box[s][MAX] = FLT_MIN;
    }
    for (int n = 0; n < poly->facetN; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            /* bounding box */
            poly->box[s][MIN] = MinReal(poly->box[s][MIN], MinReal(poly->facet[n].P1[s], MinReal(poly->facet[n].P2[s], poly->facet[n].P3[s])));
            poly->box[s][MAX] = MaxReal(poly->box[s][MAX], MaxReal(poly->facet[n].P1[s], MaxReal(poly->facet[n].P2[s], poly->facet[n].P3[s])));
            /* edge vectors */
            P1P2[s] = poly->facet[n].P2[s] - poly->facet[n].P1[s];
            P1P3[s] = poly->facet[n].P3[s] - poly->facet[n].P1[s];
        }
        /* normal vector */
        Cross(poly->facet[n].N, P1P2, P1P3);
        poly->facet[n].s = 0.5 * Norm(poly->facet[n].N);
        /* accumulate area and volume */
        poly->area = poly->area + poly->facet[n].s;
        poly->volume = poly->volume + Dot(poly->facet[n].P1, poly->facet[n].N);
        /* unit normal */
        Normalize(poly->facet[n].N, DIMS, 2.0 * poly->facet[n].s);
    }
    poly->volume = poly->volume / 6.0; /* final volume of the polyhedron */
    if (COLLAPSEN != space->collapsed) { /* space dimension collapsed */
        poly->area = poly->area - 2.0 * poly->volume; /* change to side area of a unit thickness polygon */
    }
    /* bounding sphere */
    for (int s = 0; s < DIMS; ++s) {
        poly->O[s] = 0.5 * (poly->box[s][MIN] + poly->box[s][MAX]);
        P1P2[s] = poly->box[s][MIN];
        P1P3[s] = poly->box[s][MAX];
    }
    poly->r = 0.5 * Dist(P1P2, P1P3);
    return 0;
}
int PointInPolyhedron(const int k, const int j, const int i, const Polyhedron *poly, const Space *space)
{
    const RealVector Pc = {PointSpace(i, X, space), PointSpace(j, Y, space), PointSpace(k, Z, space)};
    if (0 == poly->facetN) {
        if (0 > (Dist2(Pc, poly->O) - poly->r * poly->r)) {
            return 0;
        }
        return 1;
    }
    /*
     * Solve point-in-polyhedron problem for triangulated polyhedrons.
     * 
     * Angle summation test is used instead of ray crossing test due
     * to the simplicity and robustness, although more expensive.
     *
     * Carvalho, P. C. P., Cavalcanti, P. R. (1995). Point in polyhedron
     * testing using spherical polygons. Graphics Gems V, 42-49.
     * Note: the correct formula for the area of a spherical triangle is
     *          Delta = R * R * [(A + B + C) - pi]
     */
    const Real pi = acos(-1);
    RealVector PcP1 = {0.0};
    RealVector PcP2 = {0.0};
    RealVector PcP3 = {0.0};
    RealVector N1 = {0.0};
    RealVector N2 = {0.0};
    RealVector N3 = {0.0};
    Real normN1 = 0.0;
    Real normN2 = 0.0;
    Real normN3 = 0.0;
    Real sign = 0.0;
    Real area = 0.0;
    for (int n = 0; n < poly->facetN; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            PcP1[s] = poly->facet[n].P1[s] - Pc[s];
            PcP2[s] = poly->facet[n].P2[s] - Pc[s];
            PcP3[s] = poly->facet[n].P3[s] - Pc[s];
        }
        Cross(N1, PcP1, PcP2);
        Cross(N2, PcP2, PcP3);
        Cross(N3, PcP3, PcP1);
        normN1 = Norm(N1);
        normN2 = Norm(N2);
        normN3 = Norm(N3);
        sign = Sign(Dot(PcP1, poly->facet[n].N));
        area = area + sign * (2 * pi - acos(Dot(N1, N2) / (normN1 * normN2)) - acos(Dot(N2, N3) / (normN2 * normN3)) - acos(Dot(N3, N1) / (normN3 * normN1)));
    }
    if ((2 * pi < area) || (-2 * pi > area)) {
        return 0;
    }
    return 1;
}
/* a good practice: end file with a newline */

