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
#ifndef ARTRACFD_INITIALIZATION_H_ /* if this is the first definition */
#define ARTRACFD_INITIALIZATION_H_ /* a unique marker for this header file */
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
 * Flow field initializer
 *
 * Function
 *      Initialize the flow field variable.
 *      non restart -- initialize by specified data.
 *      restart -- initialize by restart data.
 */
extern int InitializeFlowField(Real *U, const Space *, const Geometry *, Time *,
        const Partition *, const Flow *);
#endif
/* a good practice: end file with a newline */

