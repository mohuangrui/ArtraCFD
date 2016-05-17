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
static int FindEdge(const int, const int, const Polyhedron *);
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
static int FindEdge(const int v0, const int v1, const Polyhedron *poly)
{
    /* search the edge list and return the edge index */
    for (int n = 0; n < poly->edgeN; ++n) {
        if (((v0 == poly->e[n][0]) && (v1 == poly->e[n][1])) ||
                ((v1 == poly->e[n][0]) && (v0 == poly->e[n][1]))) {
            return n;
        }
    }
    /* otherwise, something is wrong */
    return 0;
}
void Transformation(const Real O[restrict], const Real scale[restrict], const Real angle[restrict],
        const Real offset[restrict], Polyhedron *poly)
{
    const RealVec Sin = {sin(angle[X]), sin(angle[Y]), sin(angle[Z])};
    const RealVec Cos = {cos(angle[X]), cos(angle[Y]), cos(angle[Z])};
    const Real rotate[DIMS][DIMS] = {
        {Cos[Y]*Cos[Z], Cos[X]*Sin[Z]+Sin[X]*Sin[Y]*Cos[Z], Sin[X]*Sin[Z]-Cos[X]*Sin[Y]*Cos[Z]},
        {-Cos[Y]*Sin[Z], Cos[X]*Cos[Z]-Sin[X]*Sin[Y]*Sin[Z], Sin[X]*Cos[Z]+Cos[X]*Sin[Y]*Sin[Z]},
        {Sin[Y], -Sin[X]*Cos[Y], Cos[X]*Cos[Y]}};
    RealVec tmp = {0.0};
    for (int n = 0; n < poly->vertN; ++n) {
        /* move to relative origin */
        poly->v[n][X] = poly->v[n][X] - O[X];
        poly->v[n][Y] = poly->v[n][Y] - O[Y];
        poly->v[n][Z] = poly->v[n][Z] - O[Z];
        /* scale */
        poly->v[n][X] = poly->v[n][X] * scale[X];
        poly->v[n][Y] = poly->v[n][Y] * scale[Y];
        poly->v[n][Z] = poly->v[n][Z] * scale[Z];
        /* rotate */
        tmp[X] = Dot(rotate[X], poly->v[n]);
        tmp[Y] = Dot(rotate[Y], poly->v[n]);
        tmp[Z] = Dot(rotate[Z], poly->v[n]);
        /* translate */
        poly->v[n][X] = tmp[X] + offset[X] + O[X];
        poly->v[n][Y] = tmp[Y] + offset[Y] + O[Y];
        poly->v[n][Z] = tmp[Z] + offset[Z] + O[Z];
    }
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
        BuildTriangle(n, poly, v0, v1, v2, e01, e02);
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
        const RealVec e12 = {v2[X] - v1[X], v2[Y] - v1[Y], v2[Z] - v1[Z]};
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
void BuildTriangle(const int faceID, const Polyhedron *poly, Real v0[restrict], 
        Real v1[restrict], Real v2[restrict], Real e01[restrict], Real e02[restrict])
{
    for (int s = 0; s < DIMS; ++s) {
        /* vertices */
        v0[s] = poly->v[poly->f[faceID][0]][s];
        v1[s] = poly->v[poly->f[faceID][1]][s];
        v2[s] = poly->v[poly->f[faceID][2]][s];
        /* edge vectors */
        e01[s] = v1[s] - v0[s];
        e02[s] = v2[s] - v0[s];
    }
    return;
}
int PointInPolyhedron(const Real p[restrict], const Polyhedron *poly, int faceID[restrict])
{
    RealVec v0 = {0.0}; /* vertices */
    RealVec v1 = {0.0};
    RealVec v2 = {0.0};
    RealVec e01 = {0.0}; /* edges */
    RealVec e02 = {0.0};
    RealVec pi = {0.0}; /* closest point */
    RealVec N = {0.0}; /* normal of the closest point */
    /*
     * Parametric equation of triangle defined plane 
     * T(s,t) = v0 + s(v1-v0) + t(v2-v0) = v0 + s*e01 + t*e02
     * s, t: real numbers. v0, v1, v2: vertices. e01, e02: edge vectors. 
     * A point pi = T(s,t) is in the triangle T when s>=0, t>=0, and s+t<=1.
     * Further, pi is on an edge of T if one of the conditions s=0; t=0; 
     * s+t=1 is true with each condition corresponds to one edge. Each 
     * s=0, t=0; s=1, t=0; s=0, t=1 corresponds to v0, v1, and v2.
     */
    RealVec para = {0.0}; /* parametric coordinates */
    Real distSquare = 0.0; /* store computed squared distance */
    Real distSquareMin = FLT_MAX; /* store minimum squared distance */
    int fid = 0; /* closest face identifier */
    for (int n = 0; n < poly->faceN; ++n) {
        BuildTriangle(n, poly, v0, v1, v2, e01, e02);
        distSquare = PointTriangleDistance(p, v0, e01, e02, para);
        if (distSquareMin > distSquare) {
            distSquareMin = distSquare;
            fid = n;
        }
    }
    *faceID = fid;
    ComputeIntersection(p, fid, poly, pi, N);
    pi[X] = p[X] - pi[X];
    pi[Y] = p[Y] - pi[Y];
    pi[Z] = p[Z] - pi[Z];
    if (0.0 < Dot(pi, N)) {
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
    Real distSquare = 0.0;
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
Real ComputeIntersection(const Real p[restrict], const int faceID, const Polyhedron *poly,
        Real pi[restrict], Real N[restrict])
{
    RealVec v0 = {0.0}; /* vertices */
    RealVec v1 = {0.0};
    RealVec v2 = {0.0};
    RealVec e01 = {0.0}; /* edges */
    RealVec e02 = {0.0};
    RealVec para = {0.0}; /* parametric coordinates */
    const Real zero = 0.0;
    const Real one = 1.0;
    const IntVec v = {poly->f[faceID][0], poly->f[faceID][1], poly->f[faceID][2]}; /* vertex index in vertex list */
    int e = 0; /* edge index in edge list */
    BuildTriangle(faceID, poly, v0, v1, v2, e01, e02);
    const Real distSquare = PointTriangleDistance(p, v0, e01, e02, para);
    if (EqualReal(para[1], zero)) {
        if (EqualReal(para[2], zero)) {
            /* vertex 0 */
            for (int s = 0; s < DIMS; ++s) {
                pi[s] = v0[s];
                N[s] = poly->Nv[v[0]][s];
            }
        } else {
            if (EqualReal(para[2], one)) {
                /* vertex 2 */
                for (int s = 0; s < DIMS; ++s) {
                    pi[s] = v2[s];
                    N[s] = poly->Nv[v[2]][s];
                }
            } else {
                /* edge e02 */
                e = FindEdge(v[0], v[2], poly);
                for (int s = 0; s < DIMS; ++s) {
                    pi[s] = v0[s] + para[2] * e02[s];
                    N[s] = poly->Ne[e][s];
                }
            }
        }
    } else {
        if (EqualReal(para[1], one)) {
            /* vertex 1 */
            for (int s = 0; s < DIMS; ++s) {
                pi[s] = v1[s];
                N[s] = poly->Nv[v[1]][s];
            }
        } else {
            if (EqualReal(para[2], zero)) {
                /* edge e01 */
                e = FindEdge(v[0], v[1], poly);
                for (int s = 0; s < DIMS; ++s) {
                    pi[s] = v0[s] + para[1] * e01[s];
                    N[s] = poly->Ne[e][s];
                }
            } else {
                if (EqualReal(para[0], zero)) {
                    /* edge e12 */
                    e = FindEdge(v[1], v[2], poly);
                    for (int s = 0; s < DIMS; ++s) {
                        pi[s] = v0[s] + para[1] * e01[s] + para[2] * e02[s];
                        N[s] = poly->Ne[e][s];
                    }
                } else {
                    /* complete in the triangle */
                    for (int s = 0; s < DIMS; ++s) {
                        pi[s] = v0[s] + para[1] * e01[s] + para[2] * e02[s];
                        N[s] = poly->Nf[faceID][s];
                    }
                }
            }
        }
    }
    return distSquare;
}
/* a good practice: end file with a newline */

