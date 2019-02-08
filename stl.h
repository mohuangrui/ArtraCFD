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
#ifndef ARTRACFD_STL_H_ /* if undefined */
#define ARTRACFD_STL_H_ /* set a unique marker */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * STL Reader
 */
extern void ReadStlFile(const char *fname, Polyhedron *);
/*
 * STL Writer
 */
extern void WriteStlFile(const char *fname, const Polyhedron *);
#endif
/* a good practice: end file with a newline */

