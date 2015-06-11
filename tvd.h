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
#ifndef ARTRACFD_TVD_H_ /* if this is the first definition */
#define ARTRACFD_TVD_H_ /* a unique marker for this header file */
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
 * Spatial discretization and computation
 *
 * Function
 *      Compute field data to next time step based on the inputed data.
 * Parameters
 *      U -- current data, and store updated field data after computation.
 *      dt -- time step size
 *      Uswap -- a auxiliary storage space of the same size of U.
 */
extern int SpatialDiscretizationAndComputation(Real *U, const Real dt, Real *Uswap, 
        const Space *, const Geometry *, const Partition *, const Flow *);
#endif
/* a good practice: end file with a newline */

