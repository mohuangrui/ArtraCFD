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
#ifndef ARTRACFD_WENO_H_ /* if this is the first definition */
#define ARTRACFD_WENO_H_ /* a unique marker for this header file */
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
 * WENO
 *
 * Function
 *      Compute reconstructed convective fluxes by WENO scheme.
 */
extern void WENO3(const int tn, const int s, const int k, const int j, const int i, 
        const int partn[restrict], const Node *const, const Model *, Real Fhat[restrict]);
extern void WENO5(const int tn, const int s, const int k, const int j, const int i, 
        const int partn[restrict], const Node *const, const Model *, Real Fhat[restrict]);
#endif
/* a good practice: end file with a newline */

