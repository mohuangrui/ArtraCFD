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
 * Header File Guards to Avoid Interdependence
 ****************************************************************************/
#ifndef ARTRACFD_COMPUTATIONAL_GEOMETRY_H_ /* if undefined */
#define ARTRACFD_COMPUTATIONAL_GEOMETRY_H_ /* set a unique marker */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Polyhedron representation
 *
 * Function
 *      Convert polyhedron representation from STL to a mixture form
 *      of face-vertex mesh and winged-edge mesh.
 */
extern void ConvertPolyhedron(Polyhedron *);
extern void AllocatePolyhedronMemory(const int vertN, const int edgeN,
        const int faceN, Polyhedron *);
extern void AddEdge(const int v0, const int v1, const int f, Polyhedron *);
extern void QuickSortEdge(const int n, int e[restrict][EVF]);
extern void BuildTriangle(const int fid, const Polyhedron *, Real v0[restrict],
        Real v1[restrict], Real v2[restrict], Real e01[restrict], Real e02[restrict]);
/*
 * Compute geometry parameters
 *
 * Function
 *      Compute the geometric properties of each polyhedron, including bounding
 *      volume, area, volume, centroid, inertia tensor, normal. Note that the
 *      inertia tensor is relative to the body coordinates located at centroid
 *      and is computed by assuming that the density is a constant with value 1.
 */
extern void ComputeGeometryParameters(const int collapse, Geometry *const);
/*
 * Polyhedron transformation
 */
extern void TransformPolyhedron(const Real O[restrict], const Real scale[restrict],
        const Real angle[restrict], const Real offset[restrict], Polyhedron *);
/*
 * Point in polyhedron
 *
 * Function
 *      Solve point-in-polyhedron problem for triangulated polyhedron,
 *      also find the cloest face.
 */
extern int PointInPolyhedron(const Real p[restrict], const Polyhedron *, int fid[restrict]);
/*
 * Point triangle distance
 *
 * Function
 *     Returns the squared minimum distance from a point to a triangle,
 *     also finds the barycentric coordnates of the intersection point.
 */
extern Real PointTriangleDistance(const Real p[restrict], const Real v0[restrict],
        const Real e01[restrict], const Real e02[restrict], Real para[restrict]);
/*
 * Point triangle intersection point
 *
 * Function
 *      Obtain the coordinates and normal of the intersection point,
 *      also return the distance.
 */
extern Real ComputeIntersection(const Real p[restrict], const int fid,
        const Polyhedron *poly, Real pi[restrict], Real N[restrict]);
/*
 * Compute geometric data
 *
 * Function
 *      Compute the intersection point pi, mirror point pm, outward surface
 *      normal N of the point p regarding the face fid of the polyhedron.
 */
extern void ComputeGeometricData(const Real p[restrict], const int fid, const Polyhedron *,
        Real pi[restrict], Real pm[restrict], Real N[restrict]);
#endif
/* a good practice: end file with a newline */

