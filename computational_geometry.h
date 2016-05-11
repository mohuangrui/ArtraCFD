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
#ifndef ARTRACFD_COMPUTATIONAL_GEOMETRY_H_ /* if this is the first definition */
#define ARTRACFD_COMPUTATIONAL_GEOMETRY_H_ /* a unique marker for this header file */
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
extern void AllocatePolyhedronMemory(const int vertN, const int faceN, Polyhedron *);
extern void AddEdge(const int v0, const int v1, const int f, Polyhedron *);
/*
 * Compute geometry parameters
 *
 * Function
 *      Compute the geometric properties of each polyhedron.
 */
extern void ComputeGeometryParameters(const int collapse, Geometry *);
/*
 * Point in polyhedron
 *
 * Function
 *      Solve point-in-polyhedron problem for triangulated polyhedron.
 */
extern int PointInPolyhedron(const Real pc[restrict], const Polyhedron *);
#endif
/* a good practice: end file with a newline */

