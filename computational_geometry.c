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
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <string.h> /* manipulating strings */
#include <math.h> /* common mathematical functions */
#include <float.h> /* size of floating point values */
#include "cfd_commons.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int AddVertex(const Real [restrict], Polyhedron *);
static int FindEdge(const int, const int, const int, int [restrict][EVF]);
static void ComputeParametersSphere(const int, Polyhedron *);
static void ComputeParametersPolyhedron(const int, Polyhedron *);
static void TransformVertex(const Real [restrict], const Real [restrict],
        const Real [restrict][DIMS], const Real [restrict], Real [restrict][LIMIT],
        const int, Real [restrict][DIMS]);
static void TransformNormal(const Real [restrict][DIMS], const int, Real [restrict][DIMS]);
static Real TransformInertia(const Real [restrict], Real [restrict][DIMS]);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void ConvertPolyhedron(Polyhedron *poly)
{
    /* allocate memory, assume over-estimated vertex and edge */
    AllocatePolyhedronMemory(POLYN * poly->faceN, POLYN * poly->faceN, poly->faceN, poly);
    /* convert representation */
    for (int n = 0; n < poly->faceN; ++n) {
        poly->f[n][0] = AddVertex(poly->facet[n].v0, poly);
        poly->f[n][1] = AddVertex(poly->facet[n].v1, poly);
        poly->f[n][2] = AddVertex(poly->facet[n].v2, poly);
        AddEdge(poly->f[n][0], poly->f[n][1], n, poly);
        AddEdge(poly->f[n][1], poly->f[n][2], n, poly);
        AddEdge(poly->f[n][2], poly->f[n][0], n, poly);
    }
    QuickSortEdge(poly->edgeN, poly->e);
    /* adjust the memory allocation */
    RetrieveStorage(poly->facet);
    poly->facet = NULL;
    poly->e = realloc(poly->e, poly->edgeN * sizeof(*poly->e));
    poly->Ne = realloc(poly->Ne, poly->edgeN * sizeof(*poly->Ne));
    poly->v = realloc(poly->v, poly->vertN * sizeof(*poly->v));
    poly->Nv = realloc(poly->Nv, poly->vertN * sizeof(*poly->Nv));
    return;
}
void AllocatePolyhedronMemory(const int vertN, const int edgeN,
        const int faceN, Polyhedron *poly)
{
    poly->f = AssignStorage(faceN * sizeof(*poly->f));
    poly->Nf = AssignStorage(faceN * sizeof(*poly->Nf));
    poly->e = AssignStorage(edgeN * sizeof(*poly->e));
    poly->Ne = AssignStorage(edgeN * sizeof(*poly->Ne));
    poly->v = AssignStorage(vertN * sizeof(*poly->v));
    poly->Nv = AssignStorage(vertN * sizeof(*poly->Nv));
    return;
}
static int AddVertex(const Real v[restrict], Polyhedron *poly)
{
    /* search the vertex list, if already exist, return the index */
    for (int n = 0; n < poly->vertN; ++n) {
        if ((v[X] == poly->v[n][X]) && (v[Y] == poly->v[n][Y]) &&
                (v[Z] == poly->v[n][Z])) {
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
void AddEdge(const int v0, const int v1, const int f, Polyhedron *poly)
{
    /* insert by a predefined order */
    const int vMax = (v0 > v1) ? v0 : v1;
    const int vMin = (v0 > v1) ? v1 : v0;
    /* search the edge list, if already exist, add the second face index */
    for (int n = 0; n < poly->edgeN; ++n) {
        if ((vMax == poly->e[n][0]) && (vMin == poly->e[n][1])) {
            poly->e[n][3] = f;
            return;
        }
    }
    /* otherwise, add to the edge list */
    poly->e[poly->edgeN][0] = vMax;
    poly->e[poly->edgeN][1] = vMin;
    poly->e[poly->edgeN][2] = f;
    ++(poly->edgeN); /* increase pointer */
    return;
}
void QuickSortEdge(const int n, int e[restrict][EVF])
{
    if (2 > n) {
        return;
    }
    int v[2] = {0};
    int temp = 0;
    int i = 0;
    int j = 0;
    v[0] = e[n/2][0];
    v[1] = e[n/2][1];
    for (i = 0, j = n - 1;; ++i, --j) {
        while ((e[i][0] < v[0]) || ((e[i][0] == v[0]) && (e[i][1] < v[1]))) {
            ++i;
        }
        while ((e[j][0] > v[0]) || ((e[j][0] == v[0]) && (e[j][1] > v[1]))) {
            --j;
        }
        if (i >= j) {
            break;
        }
        for (int k = 0; k < EVF; ++k) {
            temp = e[i][k];
            e[i][k] = e[j][k];
            e[j][k] = temp;
        }
    }
    QuickSortEdge(i, e);
    QuickSortEdge(n - i, e + i);
    return;
}
static int FindEdge(const int v0, const int v1, const int n, int e[restrict][EVF])
{
    /* obtain a predefined order */
    const int vMax = (v0 > v1) ? v0 : v1;
    const int vMin = (v0 > v1) ? v1 : v0;
    /* binary search the edge list and return the edge index */
    int i = 0;
    int j = n - 1;
    int k = 0;
    while (i <= j) {
        k = (i + j) / 2;
        if ((vMax == e[k][0]) && (vMin == e[k][1])) {
            return k;
        } else {
            if ((vMax > e[k][0]) || ((vMax == e[k][0]) && (vMin > e[k][1]))) {
                i = k + 1;
            } else {
                j = k - 1;
            }
        }
    }
    /* target was not found */
    ShowError("finding edge failed...");
    return -1;
}
void TransformPolyhedron(const Real O[restrict], const Real scale[restrict],
        const Real angle[restrict], const Real offset[restrict], Polyhedron *poly)
{
    const RealVec Sin = {sin(angle[X]), sin(angle[Y]), sin(angle[Z])};
    const RealVec Cos = {cos(angle[X]), cos(angle[Y]), cos(angle[Z])};
    const Real rotate[DIMS][DIMS] = { /* point rotation matrix */
        {Cos[Y]*Cos[Z], -Cos[X]*Sin[Z]+Sin[X]*Sin[Y]*Cos[Z], Sin[X]*Sin[Z]+Cos[X]*Sin[Y]*Cos[Z]},
        {Cos[Y]*Sin[Z], Cos[X]*Cos[Z]+Sin[X]*Sin[Y]*Sin[Z], -Sin[X]*Cos[Z]+Cos[X]*Sin[Y]*Sin[Z]},
        {-Sin[Y], Sin[X]*Cos[Y], Cos[X]*Cos[Y]}};
    const Real invrot[DIMS][DIMS] = { /* inverse rotation matrix */
        {rotate[0][0], rotate[1][0], rotate[2][0]},
        {rotate[0][1], rotate[1][1], rotate[2][1]},
        {rotate[0][2], rotate[1][2], rotate[2][2]}};
    const Real num = 1.0 / sqrt(2.0);
    const Real axe[6][DIMS] = { /* direction vector of axis xx, yy, zz, xy, yz, zx */
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
        {num, num, 0.0}, {0.0, num, num}, {num, 0.0, num}};
    RealVec axis = {0.0}; /* direction vector of axis in rotated frame */
    Real I[6] = {0.0}; /* inertia tensor after rotation */
    /* transforming vertex and build the new bounding box */
    for (int s = 0; s < DIMS; ++s) {
        poly->box[s][MIN] = FLT_MAX;
        poly->box[s][MAX] = FLT_MIN;
    }
    TransformVertex(O, scale, rotate, offset, poly->box, poly->vertN, poly->v);
    /* transforming normal assuming pure rotation and translation */
    TransformNormal(rotate, poly->faceN, poly->Nf);
    TransformNormal(rotate, poly->edgeN, poly->Ne);
    TransformNormal(rotate, poly->vertN, poly->Nv);
    /* transform inertial tensor */
    for (int n = 0; n < 6; ++n) {
        axis[X] = Dot(invrot[X], axe[n]);
        axis[Y] = Dot(invrot[Y], axe[n]);
        axis[Z] = Dot(invrot[Z], axe[n]);
        I[n] = TransformInertia(axis, poly->I);
    }
    poly->I[X][X] = I[0];  poly->I[X][Y] = -I[3]; poly->I[X][Z] = -I[5];
    poly->I[Y][X] = -I[3]; poly->I[Y][Y] = I[1];  poly->I[Y][Z] = -I[4];
    poly->I[Z][X] = -I[5]; poly->I[Z][Y] = -I[4]; poly->I[Z][Z] = I[2];
    /* centroid should be transformed at last */
    Real Oc[1][DIMS] = {{poly->O[X], poly->O[Y], poly->O[Z]}};
    TransformVertex(O, scale, rotate, offset, poly->box, 1, Oc);
    poly->O[X] = Oc[0][X];
    poly->O[Y] = Oc[0][Y];
    poly->O[Z] = Oc[0][Z];
    return;
}
static void TransformVertex(const Real O[restrict], const Real scale[restrict],
        const Real rotate[restrict][DIMS], const Real offset[restrict],
        Real box[restrict][LIMIT], const int vertN, Real v[restrict][DIMS])
{
    RealVec tmp = {0.0};
    for (int n = 0; n < vertN; ++n) {
        /* translate reference frame to a parallel frame at the reference point */
        v[n][X] = v[n][X] - O[X];
        v[n][Y] = v[n][Y] - O[Y];
        v[n][Z] = v[n][Z] - O[Z];
        /* scale */
        v[n][X] = v[n][X] * scale[X];
        v[n][Y] = v[n][Y] * scale[Y];
        v[n][Z] = v[n][Z] * scale[Z];
        /* rotate */
        tmp[X] = v[n][X];
        tmp[Y] = v[n][Y];
        tmp[Z] = v[n][Z];
        v[n][X] = Dot(rotate[X], tmp);
        v[n][Y] = Dot(rotate[Y], tmp);
        v[n][Z] = Dot(rotate[Z], tmp);
        /* translate with offset and then translate reference back to origin */
        v[n][X] = v[n][X] + offset[X] + O[X];
        v[n][Y] = v[n][Y] + offset[Y] + O[Y];
        v[n][Z] = v[n][Z] + offset[Z] + O[Z];
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = (box[s][MIN] < v[n][s]) ? box[s][MIN] : v[n][s];
            box[s][MAX] = (box[s][MAX] > v[n][s]) ? box[s][MAX] : v[n][s];
        }
    }
    return;
}
static void TransformNormal(const Real matrix[restrict][DIMS],
        const int normalN, Real N[restrict][DIMS])
{
    RealVec tmp = {0.0};
    for (int n = 0; n < normalN; ++n) {
        tmp[X] = N[n][X];
        tmp[Y] = N[n][Y];
        tmp[Z] = N[n][Z];
        N[n][X] = Dot(matrix[X], tmp);
        N[n][Y] = Dot(matrix[Y], tmp);
        N[n][Z] = Dot(matrix[Z], tmp);
        /* normalization is needed if anisotropic transformation happens */
    }
    return;
}
static Real TransformInertia(const Real axis[restrict], Real I[restrict][DIMS])
{
    return I[X][X] * axis[X] * axis[X] + I[Y][Y] * axis[Y] * axis[Y] +
        I[Z][Z] * axis[Z] * axis[Z] + 2.0 * I[X][Y] * axis[X] * axis[Y] +
        2.0 * I[Y][Z] * axis[Y] * axis[Z] + 2.0 * I[Z][X] * axis[Z] * axis[X];
}
void ComputeGeometryParameters(const int collapse, Geometry *const geo)
{
    for (int n = 0; n < geo->sphN; ++n) {
        ComputeParametersSphere(collapse, geo->poly + n);
    }
    for (int n = geo->sphN; n < geo->totN; ++n) {
        ComputeParametersPolyhedron(collapse, geo->poly + n);
    }
    return;
}
/*
 * A bounding box and a bounding sphere are both used as bounding containers
 * to enclose a finite geometric object. Meanwhile, triangulated polyhedrons
 * and analytical spheres are unified by the using of bounding container,
 * since an analytical sphere is the bounding sphere of itself. Moreover,
 * a polyhedron with a unit length thickness is used to represent a polygon
 * with the same cross-section shape.
 */
static void ComputeParametersSphere(const int collapse, Polyhedron *poly)
{
    const Real pi = PI;
    /* bounding box */
    for (int s = 0; s < DIMS; ++s) {
        poly->box[s][MIN] = poly->O[s] - poly->r;
        poly->box[s][MAX] = poly->O[s] + poly->r;
    }
    /* geometric property */
    Real num = 0.0;
    if (COLLAPSEN == collapse) { /* no space dimension collapsed */
        poly->area = 4.0 * pi * poly->r * poly->r; /* area of a sphere */
        poly->volume = (4.0 / 3.0) * pi * poly->r * poly->r * poly->r; /* volume of a sphere */
        num = 0.4;
    } else {
        poly->area = 2.0 * pi * poly->r; /* side area of a unit thickness cylinder */
        poly->volume = pi * poly->r * poly->r; /* volume of a unit thickness cylinder */
        num = 0.5;
    }
    num = num * poly->r * poly->r * poly->volume;
    poly->I[X][X] = num;  poly->I[X][Y] = 0.0;  poly->I[X][Z] = 0.0;
    poly->I[Y][X] = 0.0;  poly->I[Y][Y] = num;  poly->I[Y][Z] = 0.0;
    poly->I[Z][X] = 0.0;  poly->I[Z][Y] = 0.0;  poly->I[Z][Z] = num;
    return;
}
static void ComputeParametersPolyhedron(const int collapse, Polyhedron *poly)
{
    /* initialize parameters */
    const Real pi = PI;
    RealVec v0 = {0.0}; /* vertices */
    RealVec v1 = {0.0};
    RealVec v2 = {0.0};
    RealVec e01 = {0.0}; /* edges */
    RealVec e02 = {0.0};
    RealVec Nf = {0.0}; /* outward normal */
    RealVec Angle = {0.0}; /* internal angle */
    RealVec O = {0.0}; /* centroid */
    Real area = 0.0; /* area */
    Real volume = 0.0; /* volume */
    Real I[6] = {0.0}; /* inertia integration xx, yy, zz, xy, yz, zx */
    RealVec tmp = {0.0}; /* temporary */
    Real f[DIMS][DIMS] = {{0.0}}; /* temporary */
    Real g[DIMS][DIMS] = {{0.0}}; /* temporary */
    Real box[LIMIT][DIMS] = {{0.0}}; /* bounding box */
    for (int s = 0; s < DIMS; ++s) {
        box[MIN][s] = FLT_MAX;
        box[MAX][s] = FLT_MIN;
    }
    /* initialize vertices normal */
    memset(poly->Nv, 0, poly->vertN * sizeof(*poly->Nv));
    /* bounding box */
    for (int n = 0; n < poly->vertN; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            box[MIN][s] = (box[MIN][s] < poly->v[n][s]) ? box[MIN][s] : poly->v[n][s];
            box[MAX][s] = (box[MAX][s] > poly->v[n][s]) ? box[MAX][s] : poly->v[n][s];
        }
    }
    /*
     * Gelder, A. V. (1995). Efficient computation of polygon area and
     * polyhedron volume. Graphics Gems V.
     * Eberly, David. "Polyhedral mass properties (revisited)." Geometric
     * Tools, LLC, Tech. Rep (2002).
     */
    for (int n = 0; n < poly->faceN; ++n) {
        BuildTriangle(n, poly, v0, v1, v2, e01, e02);
        /* outward normal vector */
        Cross(e01, e02, Nf);
        /* temporary values */
        for (int s = 0; s < DIMS; ++s) {
            tmp[0] = v0[s] + v1[s];
            f[0][s] = tmp[0] + v2[s];
            tmp[1] = v0[s] * v0[s];
            tmp[2] = tmp[1] + v1[s] * tmp[0];
            f[1][s] = tmp[2] + v2[s] * f[0][s];
            f[2][s] = v0[s] * tmp[1] + v1[s] * tmp[2] + v2[s] * f[1][s];
            g[0][s] = f[1][s] + v0[s] * (f[0][s] + v0[s]);
            g[1][s] = f[1][s] + v1[s] * (f[0][s] + v1[s]);
            g[2][s] = f[1][s] + v2[s] * (f[0][s] + v2[s]);
        }
        /* integration */
        area = area + Norm(Nf);
        volume = volume + Nf[X] * f[0][X];
        O[X] = O[X] + Nf[X] * f[1][X];
        O[Y] = O[Y] + Nf[Y] * f[1][Y];
        O[Z] = O[Z] + Nf[Z] * f[1][Z];
        I[0] = I[0] + Nf[X] * f[2][X];
        I[1] = I[1] + Nf[Y] * f[2][Y];
        I[2] = I[2] + Nf[Z] * f[2][Z];
        I[3] = I[3] + Nf[X] * (v0[Y] * g[0][X] + v1[Y] * g[1][X] + v2[Y] * g[2][X]);
        I[4] = I[4] + Nf[Y] * (v0[Z] * g[0][Y] + v1[Z] * g[1][Y] + v2[Z] * g[2][Y]);
        I[5] = I[5] + Nf[Z] * (v0[X] * g[0][Z] + v1[X] * g[1][Z] + v2[X] * g[2][Z]);
        /* unit normal */
        Normalize(DIMS, Norm(Nf), Nf);
        /*
         * Refine vertices normal by corresponding angles
         * Baerentzen, J. A., & Aanaes, H. (2005). Signed distance computation
         * using the angle weighted pseudonormal. Visualization and Computer
         * Graphics, IEEE Transactions on, 11(3), 243-253.
         */
        /* calculate internal angles by the law of cosines */
        const RealVec e12 = {v2[X] - v1[X], v2[Y] - v1[Y], v2[Z] - v1[Z]};
        const RealVec lsq = {Dot(e01, e01), Dot(e02, e02), Dot(e12, e12)};
        Angle[0] = acos((lsq[0] + lsq[1] - lsq[2]) / (2.0 * sqrt(lsq[0] * lsq[1])));
        Angle[1] = acos((lsq[0] + lsq[2] - lsq[1]) / (2.0 * sqrt(lsq[0] * lsq[2])));
        Angle[2] = pi - Angle[0] - Angle[1];
        for (int v = 0; v < POLYN; ++v) {
            for (int s = 0; s < DIMS; ++s) {
                poly->Nv[poly->f[n][v]][s] = poly->Nv[poly->f[n][v]][s] + Angle[v] * Nf[s];
            }
        }
        /* assign face normal */
        for (int s = 0; s < DIMS; ++s) {
            poly->Nf[n][s] = Nf[s];
        }
    }
    /* rectify final integration */
    area = area * (1.0 / 2.0);
    volume = volume * (1.0 / 6.0);
    O[X] = O[X] * (1.0 / 24.0);
    O[Y] = O[Y] * (1.0 / 24.0);
    O[Z] = O[Z] * (1.0 / 24.0);
    I[0] = I[0] * (1.0 / 60.0);
    I[1] = I[1] * (1.0 / 60.0);
    I[2] = I[2] * (1.0 / 60.0);
    I[3] = I[3] * (1.0 / 120.0);
    I[4] = I[4] * (1.0 / 120.0);
    I[5] = I[5] * (1.0 / 120.0);
    O[X] = O[X] / volume;
    O[Y] = O[Y] / volume;
    O[Z] = O[Z] / volume;
    /* assign to polyhedron */
    if (COLLAPSEN == collapse) { /* no space dimension collapsed */
        poly->area = area;
    } else {
        poly->area = area - 2.0 * volume; /* change to side area of a unit thickness polygon */
    }
    poly->volume = volume;
    poly->O[X] = O[X];
    poly->O[Y] = O[Y];
    poly->O[Z] = O[Z];
    /* inertia relative to centroid */
    poly->I[X][X] = I[1] + I[2] - volume * (O[Y] * O[Y] + O[Z] * O[Z]);
    poly->I[X][Y] = -I[3] + volume * O[X] * O[Y];
    poly->I[X][Z] = -I[5] + volume * O[Z] * O[X];
    poly->I[Y][X] = poly->I[X][Y];
    poly->I[Y][Y] = I[0] + I[2] - volume * (O[Z] * O[Z] + O[X] * O[X]);
    poly->I[Y][Z] = -I[4] + volume * O[Y] * O[Z];
    poly->I[Z][X] = poly->I[X][Z];
    poly->I[Z][Y] = poly->I[Y][Z];
    poly->I[Z][Z] = I[0] + I[1] - volume * (O[X] * O[X] + O[Y] * O[Y]);
    for (int s = 0; s < DIMS; ++s) {
        poly->box[s][MIN] = box[MIN][s];
        poly->box[s][MAX] = box[MAX][s];
    }
    /* a radius for estimating maximum velocity */
    poly->r = Dist(box[MIN], box[MAX]);
    /* normalize vertices normal */
    for (int n = 0; n < poly->vertN; ++n) {
        Normalize(DIMS, Norm(poly->Nv[n]), poly->Nv[n]);
    }
    /* compute edge normal */
    for (int n = 0; n < poly->edgeN; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            poly->Ne[n][s] = poly->Nf[poly->e[n][2]][s] + poly->Nf[poly->e[n][3]][s];
        }
        Normalize(DIMS, Norm(poly->Ne[n]), poly->Ne[n]);
    }
    return;
}
void BuildTriangle(const int fid, const Polyhedron *poly, Real v0[restrict],
        Real v1[restrict], Real v2[restrict], Real e01[restrict], Real e02[restrict])
{
    for (int s = 0; s < DIMS; ++s) {
        /* vertices */
        v0[s] = poly->v[poly->f[fid][0]][s];
        v1[s] = poly->v[poly->f[fid][1]][s];
        v2[s] = poly->v[poly->f[fid][2]][s];
        /* edge vectors */
        e01[s] = v1[s] - v0[s];
        e02[s] = v2[s] - v0[s];
    }
    return;
}
int PointInPolyhedron(const Real p[restrict], const Polyhedron *poly, int fid[restrict])
{
    const Real zero = 0.0;
    RealVec v0 = {zero}; /* vertices */
    RealVec v1 = {zero};
    RealVec v2 = {zero};
    RealVec e01 = {zero}; /* edges */
    RealVec e02 = {zero};
    RealVec pi = {zero}; /* closest point */
    RealVec N = {zero}; /* normal of the closest point */
    /*
     * Parametric equation of triangle defined plane
     * T(s,t) = v0 + s(v1-v0) + t(v2-v0) = v0 + s*e01 + t*e02
     * s, t: real numbers. v0, v1, v2: vertices. e01, e02: edge vectors.
     * A point pi = T(s,t) is in the triangle T when s>=0, t>=0, and s+t<=1.
     * Further, pi is on an edge of T if one of the conditions s=0; t=0;
     * s+t=1 is true with each condition corresponds to one edge. Each
     * s=0, t=0; s=1, t=0; s=0, t=1 corresponds to v0, v1, and v2.
     */
    RealVec para = {zero}; /* parametric coordinates */
    Real distSquare = zero; /* store computed squared distance */
    Real distSquareMin = FLT_MAX; /* store minimum squared distance */
    int cid = 0; /* closest face identifier */
    for (int n = 0; n < poly->faceN; ++n) {
        BuildTriangle(n, poly, v0, v1, v2, e01, e02);
        distSquare = PointTriangleDistance(p, v0, e01, e02, para);
        if (distSquareMin > distSquare) {
            distSquareMin = distSquare;
            cid = n;
        }
    }
    *fid = cid;
    ComputeIntersection(p, cid, poly, pi, N);
    pi[X] = p[X] - pi[X];
    pi[Y] = p[Y] - pi[Y];
    pi[Z] = p[Z] - pi[Z];
    if (zero < Dot(pi, N)) {
        /* outside polyhedron */
        return 0;
    } else {
        /* inside or on polyhedron */
        return 1;
    }
}
/*
 * Eberly, D. (1999). Distance between point and triangle in 3D.
 * http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
 */
Real PointTriangleDistance(const Real p[restrict], const Real v0[restrict], const Real e01[restrict],
        const Real e02[restrict], Real para[restrict])
{
    const RealVec D = {v0[X] - p[X], v0[Y] - p[Y], v0[Z] - p[Z]};
    const Real a = Dot(e01, e01);
    const Real b = Dot(e01, e02);
    const Real c = Dot(e02, e02);
    const Real d = Dot(e01, D);
    const Real e = Dot(e02, D);
    const Real f = Dot(D, D);
    const Real det = a * c - b * b;
    const Real zero = 0.0;
    const Real one = 1.0;
    Real s = b * e - c * d;
    Real t = b * d - a * e;
    Real distSquare = zero;
    if (s + t <= det) {
        if (s < zero) {
            if (t < zero) {
                /* region 4 */;
                if (d < zero) {
                    t = zero;
                    if (-d >= a) {
                        s = one;
                        distSquare = a + 2.0 * d + f;
                    } else {
                        s = -d / a;
                        distSquare = d * s + f;
                    }
                } else {
                    s = zero;
                    if (e >= zero) {
                        t = zero;
                        distSquare = f;
                    } else {
                        if (-e >= c) {
                            t = one;
                            distSquare = c + 2.0 * e + f;
                        } else {
                            t = -e / c;
                            distSquare = e * t + f;
                        }
                    }
                }
            } else {
                /* region 3 */
                s = zero;
                if (e >= zero) {
                    t = zero;
                    distSquare = f;
                } else {
                    if (-e >= c) {
                        t = one;
                        distSquare = c + 2.0 * e + f;
                    } else {
                        t = -e / c;
                        distSquare = e * t + f;
                    }
                }
            }
        } else {
            if (t < zero) {
                /* region 5 */;
                t = zero;
                if (d >= zero) {
                    s = zero;
                    distSquare = f;
                } else {
                    if (-d >= a) {
                        s = one;
                        distSquare = a + 2.0 * d + f;
                    } else {
                        s = -d / a;
                        distSquare = d * s + f;
                    }
                }
            } else {
                /* region 0 */
                s = s / det;
                t = t / det;
                distSquare = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
            }
        }
    } else {
        if (s < zero) {
            /* region 2 */
            const Real tmp0 = b + d;
            const Real tmp1 = c + e;
            if (tmp1 > tmp0) {
                const Real numer = tmp1 - tmp0;
                const Real denom = a - 2.0 * b + c;
                if (numer >= denom) {
                    s = one;
                    t = zero;
                    distSquare = a + 2.0 * d + f;
                } else {
                    s = numer / denom;
                    t = one - s;
                    distSquare = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
                }
            } else {
                s = zero;
                if (tmp1 <= zero) {
                    t = one;
                    distSquare = c + 2.0 * e + f;
                } else {
                    if (e >= zero) {
                        t = zero;
                        distSquare = f;
                    } else {
                        t = -e / c;
                        distSquare = e * t + f;
                    }
                }
            }
        } else {
            if (t < zero) {
                /* region 6 */;
                const Real tmp0 = b + e;
                const Real tmp1 = a + d;
                if (tmp1 > tmp0) {
                    const Real numer = tmp1 - tmp0;
                    const Real denom = a - 2.0 * b + c;
                    if (numer >= denom) {
                        t = one;
                        s = zero;
                        distSquare = c + 2.0 * e + f;
                    } else {
                        t = numer / denom;
                        s = one - t;
                        distSquare = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
                    }
                } else {
                    t = zero;
                    if (tmp1 <= zero) {
                        s = one;
                        distSquare = a + 2.0 * d + f;
                    } else {
                        if (d >= zero) {
                            s = zero;
                            distSquare = f;
                        } else {
                            s = -d / a;
                            distSquare = d * s + f;
                        }
                    }
                }
            } else {
                /* region 1 */
                const Real numer = c + e - b - d;
                if (numer <= zero) {
                    s = zero;
                    t = one;
                    distSquare = c + 2.0 * e + f;
                } else {
                    const Real denom = a - 2.0 * b + c;
                    if (numer >= denom) {
                        s = one;
                        t = zero;
                        distSquare = a + 2.0 * d + f;
                    } else {
                        s = numer / denom;
                        t = one - s;
                        distSquare = s * (a * s + b * t + 2.0 * d) + t * (b * s + c * t + 2.0 * e) + f;
                    }
                }
            }
        }
    }
    para[0] = one - s - t;
    para[1] = s;
    para[2] = t;
    /* account for numerical round-off error */
    if (zero > distSquare) {
        distSquare = zero;
    }
    return distSquare;
}
Real ComputeIntersection(const Real p[restrict], const int fid,
        const Polyhedron *poly, Real pi[restrict], Real N[restrict])
{
    const Real zero = 0.0;
    const Real one = 1.0;
    RealVec v0 = {zero}; /* vertices */
    RealVec v1 = {zero};
    RealVec v2 = {zero};
    RealVec e01 = {zero}; /* edges */
    RealVec e02 = {zero};
    RealVec para = {zero}; /* parametric coordinates */
    const IntVec v = {poly->f[fid][0], poly->f[fid][1], poly->f[fid][2]}; /* vertex index in vertex list */
    int e = 0; /* edge index in edge list */
    BuildTriangle(fid, poly, v0, v1, v2, e01, e02);
    const Real distSquare = PointTriangleDistance(p, v0, e01, e02, para);
    if (zero == para[1]) {
        if (zero == para[2]) {
            /* vertex 0 */
            for (int s = 0; s < DIMS; ++s) {
                pi[s] = v0[s];
                N[s] = poly->Nv[v[0]][s];
            }
        } else {
            if (one == para[2]) {
                /* vertex 2 */
                for (int s = 0; s < DIMS; ++s) {
                    pi[s] = v2[s];
                    N[s] = poly->Nv[v[2]][s];
                }
            } else {
                /* edge e02 */
                e = FindEdge(v[0], v[2], poly->edgeN, poly->e);
                for (int s = 0; s < DIMS; ++s) {
                    pi[s] = v0[s] + para[2] * e02[s];
                    N[s] = poly->Ne[e][s];
                }
            }
        }
    } else {
        if (one == para[1]) {
            /* vertex 1 */
            for (int s = 0; s < DIMS; ++s) {
                pi[s] = v1[s];
                N[s] = poly->Nv[v[1]][s];
            }
        } else {
            if (zero == para[2]) {
                /* edge e01 */
                e = FindEdge(v[0], v[1], poly->edgeN, poly->e);
                for (int s = 0; s < DIMS; ++s) {
                    pi[s] = v0[s] + para[1] * e01[s];
                    N[s] = poly->Ne[e][s];
                }
            } else {
                if (zero == para[0]) {
                    /* edge e12 */
                    e = FindEdge(v[1], v[2], poly->edgeN, poly->e);
                    for (int s = 0; s < DIMS; ++s) {
                        pi[s] = v0[s] + para[1] * e01[s] + para[2] * e02[s];
                        N[s] = poly->Ne[e][s];
                    }
                } else {
                    /* complete in the triangle */
                    for (int s = 0; s < DIMS; ++s) {
                        pi[s] = v0[s] + para[1] * e01[s] + para[2] * e02[s];
                        N[s] = poly->Nf[fid][s];
                    }
                }
            }
        }
    }
    return distSquare;
}
void ComputeGeometricData(const Real p[restrict], const int fid, const Polyhedron *poly,
        Real pi[restrict], Real pm[restrict], Real N[restrict])
{
    if (0 >= poly->faceN) { /* analytical polyhedron */
        Real dist = 0.0;
        N[X] = p[X] - poly->O[X];
        N[Y] = p[Y] - poly->O[Y];
        N[Z] = p[Z] - poly->O[Z];
        dist = Norm(N);
        Normalize(DIMS, dist, N);
        dist = poly->r - dist;
        pi[X] = p[X] + dist * N[X];
        pi[Y] = p[Y] + dist * N[Y];
        pi[Z] = p[Z] + dist * N[Z];
    } else { /* triangulated polyhedron */
        ComputeIntersection(p, fid, poly, pi, N);
    }
    pm[X] = pi[X] + pi[X] - p[X];
    pm[Y] = pi[Y] + pi[Y] - p[Y];
    pm[Z] = pi[Z] + pi[Z] - p[Z];
    return;
}
/* a good practice: end file with a newline */

