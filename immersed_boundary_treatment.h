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
#ifndef ARTRACFD_GCIBM_H_ /* if this is the first definition */
#define ARTRACFD_GCIBM_H_ /* a unique marker for this header file */
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
extern int ComputeGeometryDomain(Space *, const Partition *, const Geometry *);
/*
 * Boundary treatments for ghost nodes
 *
 * Function
 *      Apply boundary conditions and treatments to ghost nodes.
 */
extern int BoundaryTreatmentsGCIBM(Real *U, const Space *, const Model *,
        const Partition *, const Geometry *);
/*
 * Reconstruct flow variables by Inverse Distance Weighting
 *
 * Function
 *      Reconstruction of the values of primitive vector Uo for a spatial
 *      point (z, y, x) based on the neighbours around node (k, j, i) in
 *      index range h and a nearby boundary point by inversed distance weighting.
 */
extern int FlowReconstruction(Real Uo[], const Real z, const Real y, const Real x,
        const int k, const int j, const int i, const int h, 
        Real UoBC[], const Real info[], const Real *geo, const Real *U,
        const Space *space, const Model *model, const Geometry *geometry);
/*
 * In geometry criteria
 *
 * Function
 *      Check node whether locates in the geometry pointed by the geometry
 *      pointer.
 *
 * Return
 *      negative -- in current geometry
 *      zero     -- on current geometry
 *      positive -- out of current geometry
 */
extern Real InGeometry(const int k, const int j, const int i, const Real *geo, const Space *);
/*
 * Calculate geometry information
 *
 * Function
 *
 *      Calculate geometry information of current node.
 */
extern int CalculateGeometryInformation(Real info[], const int k, const int j, 
        const int i, const Real *geo, const Space *);
#endif
/* a good practice: end file with a newline */

