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
#ifndef ARTRACFD_PARASIGHT_STREAM_H_ /* if this is the first definition */
#define ARTRACFD_PARASIGHT_STREAM_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * Parasight data format and type control, Parasight means the Ensight data
 * format adjusted for Paraview. This adjustment is required for the following
 * reason: the handling of complex geometry in Paraview is achieved by data 
 * filter rather than the blankID of nodes. Therefore, there is no need to
 * write changing geometry files, but need to write a data entry for filter.
 * Because of the same reason, there is no need to write each subparts, since
 * they can be distinguished from each other by the data filter.
 *
 */
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Parasight format data exporter
 *
 * Function
 *      Export conservative field data vector variable:
 *      U = [rho, rho_u, rho_v, rho_w, rho_eT]
 *      to primitive variables = [rho, u, v, w, p, T]
 *      to binary data files with Parasight data format.
 * Notice
 *      U is a linear array that stores all the values.
 *      These data are in sequential state 
 *      and can be accessed by linear index math.
 */
extern int WriteComputedDataParasight(const Real * U, const Space *, 
        const Time *, const Model *, const Partition *);
/*
 * Parasight format data loader
 *
 * Function
 *      Load computed data from output files which are written in Parasight
 *      format.
 */
extern int LoadComputedDataParasight(Real *U, const Space *, Time *, const Model *,
        const Partition *);
#endif
/* a good practice: end file with a newline */

 
