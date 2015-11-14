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
#ifndef ARTRACFD_CASE_LOADER_H_ /* if this is the first definition */
#define ARTRACFD_CASE_LOADER_H_ /* a unique marker for this header file */
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
 * Case setting data loader
 *
 * Function
 *      Load case setting data from file artracfd.case, if file does not 
 *      exist or any illegal data exists in the file, it terminates the
 *      program.
 * Outcome
 *      Will write loaded data to file artracfd.verify for verification.
 */
extern int LoadCaseSettingData(Space *, Time *, Model *, Partition *);
#endif
/* a good practice: end file with a newline */

