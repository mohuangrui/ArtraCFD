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
#ifndef ARTRACFD_SOURCE_TERM_H_ /* if undefined */
#define ARTRACFD_SOURCE_TERM_H_ /* set a unique marker */
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
 * Source term
 *
 * Function
 *      Compute the source term.
 */
void ComputePhi(const int tn, const int k, const int j, const int i,
        const int partn[restrict], const Node *const,
        const Model *, Real Phi[restrict]);
#endif
/* a good practice: end file with a newline */

