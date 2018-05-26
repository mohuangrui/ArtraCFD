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
#ifndef ARTRACFD_IMMERSED_BOUNDARY_H_ /* if this is the first definition */
#define ARTRACFD_IMMERSED_BOUNDARY_H_ /* a unique marker for this header file */
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
 * Compute geometry domain
 *
 * Function
 *      Employ node flagging algorithm to handle complex geometry that locates in
 *      the computational domain.
 */
extern void ComputeGeometryDomain(Space *, const Model *);
/*
 * Compute geometric data
 */
extern void ComputeGeometricData(const int fid, const Polyhedron *, const Real pG[restrict],
        Real pO[restrict], Real pI[restrict], Real N[restrict]);
/*
 * Immersed boundary treatments
 *
 * Function
 *      Apply boundary conditions and treatments for immersed boundaries.
 */
extern void ImmersedBoundaryTreatment(const int tn, Space *, const Model *);
extern void MethodOfImage(const Real UoI[restrict], const Real UoO[restrict], Real UoG[restrict]);
#endif
/* a good practice: end file with a newline */

