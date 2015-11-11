/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) 2014-2018 Huangrui Mo <huangrui.mo@gmail.com>              *
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
static int ComputeGeometryParameters(const Space *, const Model *, Geometry *);
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

