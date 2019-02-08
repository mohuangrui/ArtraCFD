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
#ifndef ARTRACFD_IMMERSED_BOUNDARY_H_ /* if undefined */
#define ARTRACFD_IMMERSED_BOUNDARY_H_ /* set a unique marker */
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
 * Compute geometric field
 *
 * Function
 *      Employ a multidomain node mapping algorithm to map the geometry set.
 */
extern void ComputeGeometricField(Space *, const Model *);
/*
 * Immersed boundary treatments
 *
 * Function
 *      Apply boundary conditions and treatments for immersed boundaries.
 */
extern void TreatImmersedBoundary(const int tn, Space *, const Model *);
extern void DoMethodOfImage(const Real UoI[restrict], const Real UoO[restrict], Real UoG[restrict]);
#endif
/* a good practice: end file with a newline */

