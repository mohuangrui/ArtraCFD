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
#ifndef ARTRACFD_CALCULATOR_H_ /* if undefined */
#define ARTRACFD_CALCULATOR_H_ /* set a unique marker */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef struct {
    Real t;
    Real x;
    Real y;
    Real z;
    Real ans; /* store answer */
    const Real pi;
} CalcVar; /* a set of valid variables */
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Expression calculator
 */
extern int RunCalculator(void);
/*
 * Calculate expression
 *
 * Function
 *      Calculate expressions involving a set of defined variables
 */
extern Real ComputeExpression(CalcVar *, const char *str);
#endif
/* a good practice: end file with a newline */

