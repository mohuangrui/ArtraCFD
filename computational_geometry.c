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
    /* allocate memory, assume vertN = faceN */
    AllocatePolyhedronMemory(poly->faceN, poly->faceN, poly);
    /* convert representation */
    for (int n = 0; n < poly->faceN; ++n) {
        poly->f[n][0] = AddVertex(poly->facet[n].v0, poly);
        poly->f[n][1] = AddVertex(poly->facet[n].v1, poly);
        poly->f[n][2] = AddVertex(poly->facet[n].v2, poly);
        AddEdge(poly->f[n][0], poly->f[n][1], n, poly); 
        AddEdge(poly->f[n][1], poly->f[n][2], n, poly); 
        AddEdge(poly->f[n][2], poly->f[n][0], n, poly); 
    }
    /* adjust the memory allocation */
    RetrieveStorage(poly->facet);
    poly->facet = NULL;
    poly->v = realloc(poly->v, poly->vertN * sizeof(*poly->v));
    poly->Nv = realloc(poly->Nv, poly->vertN * sizeof(*poly->Nv));
    return;
}
void AllocatePolyhedronMemory(const int vertN, const int faceN, Polyhedron *poly)
{
    /* edgeN = faceN * 3 / 2 */
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
void AddEdge(const int v0, const int v1, const int f, Polyhedron *poly)
{
    /* search the edge list, if already exist, add the second face index */
    for (int n = 0; n < poly->edgeN; ++n) {
        if (((v0 == poly->e[n][0]) && (v1 == poly->e[n][1])) ||
                ((v1 == poly->e[n][0]) && (v0 == poly->e[n][1]))) {
            poly->e[n][3] = f;
            return;
        }
    }
    /* otherwise, add to the edge list */
    poly->e[poly->edgeN][0] = v0;
    poly->e[poly->edgeN][1] = v1;
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
    /* bounding box */
    for (int s = 0; s < DIMS; ++s) {
        poly->box[s][MIN] = poly->O[s] - poly->r;
        poly->box[s][MAX] = poly->O[s] + poly->r;
    }
    /* geometric property */
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
    RealVec v0 = {0.0}; /* vertices */
    RealVec v1 = {0.0};
    RealVec v2 = {0.0};
    RealVec e01 = {0.0}; /* edges */
    RealVec e02 = {0.0};
    RealVec e12 = {0.0};
    RealVec Nf = {0.0}; /* outward normal */
    RealVec Angle = {0.0}; /* internal angle */
    RealVec O = {0.0}; /* centroid */
    Real area = 0.0; /* area */
    Real volume = 0.0; /* volume */
    Real box[DIMS][LIMIT] = {{0.0}}; /* bounding box */
    for (int s = 0; s < DIMS; ++s) {
        box[s][MIN] = FLT_MAX;
        box[s][MAX] = FLT_MIN;
    }
    /* initialize vertices normal */
    for (int n = 0; n < poly->vertN; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            poly->Nv[n][s] = 0.0;
        }
    }
    /* bounding box */
    for (int n = 0; n < poly->vertN; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            box[s][MIN] = MinReal(box[s][MIN], poly->v[n][s]);
            box[s][MAX] = MaxReal(box[s][MAX], poly->v[n][s]);
        }
    }
    /*
     * Gelder, A. V. (1995). Efficient computation of polygon area and
     * polyhedron volume. Graphics Gems V.
     * Robert Nurnberg, (2013). Calculating the area and centroid of a
     * polygon and polyhedron. Imperial College London.
     */
    for (int n = 0; n < poly->faceN; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            /* vertices */
            v0[s] = poly->v[poly->f[n][0]][s];
            v1[s] = poly->v[poly->f[n][1]][s];
            v2[s] = poly->v[poly->f[n][2]][s];
            /* edge vectors */
            e01[s] = v1[s] - v0[s];
            e02[s] = v2[s] - v0[s];
            e12[s] = v2[s] - v1[s];
        }
        /* outward normal vector */
        Cross(e01, e02, Nf);
        /* accumulate area and volume */
        area = area + Norm(Nf);
        volume = volume + Dot(v0, Nf);
        /* centroid */
        for (int s = 0; s < DIMS; ++s) {
            O[s] = O[s] + Nf[s] * (Square(v0[s] + v1[s]) + Square(v1[s] + v2[s]) + Square(v2[s] + v0[s]));
        }
        /* unit normal */
        Normalize(DIMS, Norm(Nf), Nf);
        /*
         * Refine vertices normal by corresponding angles
         * Baerentzen, J. A., & Aanaes, H. (2005). Signed distance computation
         * using the angle weighted pseudonormal. Visualization and Computer
         * Graphics, IEEE Transactions on, 11(3), 243-253.
         */
        /* calculate internal angles by the law of cosines */
        Angle[0] = acos((Dot(e01, e01) + Dot(e02, e02) - Dot(e12, e12)) / 
                (2.0 * sqrt(Dot(e01, e01) * Dot(e02, e02))));
        Angle[1] = acos((Dot(e01, e01) + Dot(e12, e12) - Dot(e02, e02)) / 
                (2.0 * sqrt(Dot(e01, e01) * Dot(e12, e12))));
        Angle[2] = 3.14159265359 - Angle[0] - Angle[1];
        for (int s = 0; s < DIMS; ++s) {
            poly->Nv[poly->f[n][0]][s] = poly->Nv[poly->f[n][0]][s] + Angle[0] * Nf[s];
            poly->Nv[poly->f[n][1]][s] = poly->Nv[poly->f[n][1]][s] + Angle[1] * Nf[s];
            poly->Nv[poly->f[n][2]][s] = poly->Nv[poly->f[n][2]][s] + Angle[2] * Nf[s];
        }
        /* assign face normal */
        for (int s = 0; s < DIMS; ++s) {
            poly->Nf[n][s] = Nf[s];
        }
    }
    area = 0.5 * area; /* final area of the polyhedron */
    volume = volume / 6.0; /* final volume of the polyhedron */
    if (COLLAPSEN == collapse) { /* no space dimension collapsed */
        poly->area = area;
    } else {
        poly->area = area - 2.0 * volume; /* change to side area of a unit thickness polygon */
    }
    poly->volume = volume;
    for (int s = 0; s < DIMS; ++s) { /* assign final centroid and bounding box */
        poly->O[s] = O[s] / (48.0 * volume);
        poly->box[s][MIN] = box[s][MIN];
        poly->box[s][MAX] = box[s][MAX];
    }
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
int PointInPolyhedron(const Real p[restrict], const Polyhedron *poly)
{
    RealVec v0 = {0.0}; /* vertices */
    RealVec v1 = {0.0};
    RealVec v2 = {0.0};
    RealVec e01 = {0.0}; /* edges */
    RealVec e02 = {0.0};
    RealVec N = {0.0}; /* normal of the closet point */
    /*
     * Parametric equation of triangle defined plane 
     * T(s,t) = v0 + s(v1-v0) + t(v2-v0) = v0 + s*e01 + s*e02
     * s, t: real numbers. v0, v1, v2: vertices. e01, e02: edge vertices. 
     * A point pi = T(s,t) is in the triangle T when s>=0, t>=0, and s+t<=1.
     * Further, pi is on an edge of T if one of the conditions s=0; t=0; 
     * s+t=1 is true with each condition corresponds to one edge. Each 
     * s=0, t=0; s=1, t=0; s=0, t=1 corresponds to v0, v1, and v2.
     */
    Real s = 0.0; /* parametric coordinates */
    Real t = 0.0; /* parametric coordinates */
    Real distSquareMin = FLT_MAX; /* store minimum squared distance */
    Real distSquare = 0.0; /* store computed squared distance */
    int faceID = 0; /* cloest face identifier */
    for (int n = 0; n < poly->faceN; ++n) {
        for (int s = 0; s < DIMS; ++s) {
            /* vertices */
            v0[s] = poly->v[poly->f[n][0]][s];
            v1[s] = poly->v[poly->f[n][1]][s];
            v2[s] = poly->v[poly->f[n][2]][s];
            /* edge vectors */
            e01[s] = v1[s] - v0[s];
            e02[s] = v2[s] - v0[s];
        }
        dist = PointTriangleDistance(p, v0, e01, e02, &s, &t);
        if (distSquareMin > dist) {
            distSquareMin = dist;
            faceID = n;
        }
    }
    return faceID;
}
/*
 * Eberly, D. (1999). Distance between point and triangle in 3D.
 * http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
 */
Real PointTriangleDistance(const Real p[restrict], const Real v0[restrict], const Real e01[restrict], 
        const Real e02[restrict], Real sPara[restrict], Real tPara[restrict])
{
    const RealVec D = {v0[X] - p[X], v0[Y] - p[Y], v0[Z] - p[Z]};
    const Real a = Dot(e01, e01);
    const Real b = Dot(e01, e02);
    const Real c = Dot(e02, e02);
    const Real d = Dot(e01, D);
    const Real e = Dot(e02, D);
    const Real f = Dot(D, D);
    const Real det = a * c - b * b;
    Real s = b * e - c * d;
    Real t = b * d - a * e;
    if (s + t <= det) {
        if (s < 0.0) {
            if (t < 0.0) {
                /* region 4 */;
                if (d < 0.0) {
                    t = 0.0;
                    if (-d >= a) {
                        s = 1.0;
                    } else {
                        s = -d / a;
                    }
                } else {
                    s = 0.0;
                    if (e >= 0.0) {
                        t = 0.0;
                    } else {
                        if (-e >= c) {
                            t = 1.0;
                        } else {
                            t = -e / c;
                        }
                    }
                }
            } else {
                /* region 3 */
                s = 0.0;
                if (e >= 0.0) {
                    t = 0.0;
                } else {
                    if (-e >= c) {
                        t = 1.0;
                    } else {
                        t = -e / c;
                    }
                }
            }        
        } else {
            if (t < 0.0) {
                /* region 5 */;
                t = 0.0;
                if (d >= 0.0) {
                    s = 0.0;
                } else {
                    if (-d >= a) {
                        s = 1.0;
                    } else {
                        s = -d / a;
                    }
                }
            } else {
                /* region 0 */
                s = s / det;
                t = t / det;
            }
        }
    } else {
        if (s < 0.0) {
            /* region 2 */
            const Real tmp0 = b + d;
            const Real tmp1 = c + e;
            if (tmp1 > tmp0) {
                const Real numer = tmp1 - tmp0;
                const Real denom = a - 2.0 * b + c;
                if (numer >= denom) {
                    s = 1.0;
                    t = 0.0;
                } else {
                    s = numer / denom;
                    t = 1.0 - s;
                }
            } else {
                s = 0.0;
                if (tmp1 <= 0.0) {
                    t = 1.0;
                } else {
                    if (e >= 0.0) {
                        t = 0.0;
                    } else {
                        t = -e / c;
                    }
                }
            }
        } else {
            if (t < 0.0) {
                /* region 6 */;
                const Real tmp0 = b + e;
                const Real tmp1 = a + d;
                if (tmp1 > tmp0) {
                    const Real numer = tmp1 - tmp0;
                    const Real denom = a - 2.0 * b + c;
                    if (numer >= denom) {
                        t = 1.0;
                        s = 0.0;
                    } else {
                        t = numer / denom;
                        s = 1.0 - t;
                    }
                } else {
                    t = 0.0;
                    if (tmp1 <= 0.0) {
                        s = 1.0;
                    } else {
                        if (d >= 0.0) {
                            s = 0.0;
                        } else {
                            s = -d / a;
                        }
                    }
                }
            } else {
                /* region 1 */
                const Real numer = c + e - b - d;
                if (numer <= 0.0) {
                    s = 0.0;
                    t = 1.0;
                } else {
                    const Real denom = a - 2.0 * b + c;
                    if (numer >= denom) {
                        s = 1.0;
                        t = 0.0;
                    } else {
                        s = numer / denom;
                        t = 1.0 - s;
                    }
                    s = (numer >= denom ? 1.0 : numer / denom);
                }
            }
        }
    }
    *sPara = s;
    *tPara = t;
}
/* a good practice: end file with a newline */

