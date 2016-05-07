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
static int AddVertex(const Real [restrict], Polyhedron *);
static void ComputeParametersSphere(const int, Polyhedron *);
static void ComputeParametersPolyhedron(const int, Polyhedron *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void ConvertPolyhedron(Polyhedron *poly)
{
    AllocatePolyhedronMemory(poly->faceN, poly->faceN, poly);
    for (int n = 0; n < poly->faceN; ++n) {
        poly->f[n][0] = AddVertex(poly->facet[n].v0, poly);
        poly->f[n][1] = AddVertex(poly->facet[n].v1, poly);
        poly->f[n][2] = AddVertex(poly->facet[n].v2, poly);
        AddEdge(poly->f[n][0], poly->f[n][1], n, poly); 
        AddEdge(poly->f[n][1], poly->f[n][2], n, poly); 
        AddEdge(poly->f[n][2], poly->f[n][0], n, poly); 
    }
    RetrieveStorage(poly->facet);
    poly->facet = NULL;
    return;
}
void AllocatePolyhedronMemory(const int vertN, const int faceN, Polyhedron *poly)
{
    /* vertN <= faceN, edgeN = faceN * 3 / 2 */
    poly->f = AssignStorage(faceN * sizeof(*poly->f));
    poly->Nf = AssignStorage(faceN * sizeof(*poly->Nf));
    poly->e = AssignStorage((int)(1.5 * faceN + 0.5) * sizeof(*poly->e));
    poly->Ne = AssignStorage((int)(1.5 * faceN + 0.5) * sizeof(*poly->Ne));
    poly->v = AssignStorage(vertN * sizeof(*poly->v));
    poly->Nv = AssignStorage(vertN * sizeof(*poly->Nv));
}
static int AddVertex(const Real v[restrict], Polyhedron *poly)
{
    /* search the vertex list, if already exist, return the index */
    for (int n = 0; n < poly->vertN; ++n) {
        if (EqualReal(v[X], poly->v[n][X]) && EqualReal(v[Y], poly->v[n][Y]) && 
                EqualReal(v[Z], poly->v[n][Z])) {
            return n;
        }
    }
    /* otherwise, add to the vertex list */
    poly->v[poly->vertN][X] = v[X];
    poly->v[poly->vertN][Y] = v[Y];
    poly->v[poly->vertN][Z] = v[Z];
    ++(poly->vertN); /* increase pointer */
    return (poly->vertN - 1); /* return index */
}
void AddEdge(const int v1, const int v2, const int f, Polyhedron *poly)
{
    /* search the edge list, if already exist, add the second face index */
    for (int n = 0; n < poly->edgeN; ++n) {
        if (((v1 == poly->e[n][0]) && (v2 == poly->e[n][1])) ||
                ((v2 == poly->e[n][0]) && (v1 == poly->e[n][1]))) {
            poly->e[n][3] = f;
            return;
        }
    }
    /* otherwise, add to the edge list */
    poly->e[poly->edgeN][0] = v1;
    poly->e[poly->edgeN][1] = v2;
    poly->e[poly->edgeN][2] = f;
    ++(poly->edgeN); /* increase pointer */
    return;
}
void ComputeGeometryParameters(const int collapse, Geometry *geo)
{
    for (int n = 0; n < geo->sphereN; ++n) {
        ComputeParametersSphere(collapse, geo->poly + n);
    }
    for (int n = geo->sphereN; n < geo->totalN; ++n) {
        ComputeParametersPolyhedron(collapse, geo->poly + n);
    }
    return;
}
static void ComputeParametersSphere(const int collapse, Polyhedron *poly)
{
    const Real pi = 3.14159265359;
    for (int s = 0; s < DIMS; ++s) {
        poly->box[s][MIN] = poly->O[s] - poly->r;
        poly->box[s][MAX] = poly->O[s] + poly->r;
    }
    if (COLLAPSEN == collapse) { /* no space dimension collapsed */
        poly->area = 4.0 * pi * poly->r * poly->r; /* area of a sphere */
        poly->volume = 4.0 * pi * poly->r * poly->r * poly->r / 3.0; /* volume of a sphere */
    } else {
        poly->area = 2.0 * pi * poly->r; /* side area of a unit thickness cylinder */
        poly->volume = pi * poly->r * poly->r; /* volume of a unit thickness cylinder */
    }
    return;
}
static void ComputeParametersPolyhedron(const int collapse, Polyhedron *poly)
{
    /* initialize parameters */
    RealVec P1P2 = {0.0};
    RealVec P1P3 = {0.0};
    poly->area = 0.0;
    poly->volume = 0.0;
    for (int s = 0; s < DIMS; ++s) {
        poly->box[s][MIN] = FLT_MAX;
        poly->box[s][MAX] = FLT_MIN;
    }
    /*
     * Gelder, A. V. (1995). Efficient computation of polygon area and
     * polyhedron volume. Graphics Gems V.
     * Robert Nurnberg, (2013). Calculating the area and centroid of a
     * polygon and polyhedron. Imperial College London.
     */
    for (int n = 0; n < poly->faceN; ++n) {
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
    const RealVec Pc = {PointSpace(i, X, space), PointSpace(j, Y, space), PointSpace(k, Z, space)};
    if (0 == poly->faceN) {
        if (0 > (Dist2(Pc, poly->O) - poly->r * poly->r)) {
            return 0;
        }
        return 1;
    }
    /*
     * Solve point-in-polyhedron problem for triangulated polyhedrons.
     *
     * Haines, E. (1994). Point in polygon strategies. Graphics gems IV.
     * Intersections of rays and triangles (3D), http://geomalgorithms.com/a06-_intersect-2.html
     * Inclusion of a point in a polygon, http://geomalgorithms.com/a03-_inclusion.html
     * 
     * Angle summation test is used instead of ray crossing test due
     * to the simplicity and robustness, although more expensive.
     * Carvalho, P. C. P., Cavalcanti, P. R. (1995). Point in polyhedron
     * testing using spherical polygons. Graphics Gems V.
     * Note: the correct formula for the area of a spherical triangle is
     *             Delta = R * R * [(A + B + C) - pi]
     * Spherical Triangle, http://mathworld.wolfram.com/SphericalTriangle.html
     */
    const Real pi = acos(-1);
    RealVec PcP1 = {0.0};
    RealVec PcP2 = {0.0};
    RealVec PcP3 = {0.0};
    RealVec N1 = {0.0};
    RealVec N2 = {0.0};
    RealVec N3 = {0.0};
    Real normN1 = 0.0;
    Real normN2 = 0.0;
    Real normN3 = 0.0;
    Real sign = 0.0;
    Real area = 0.0;
    for (int n = 0; n < poly->faceN; ++n) {
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

