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
#ifndef ARTRACFD_STL_H_ /* if this is the first definition */
#define ARTRACFD_STL_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * STL data format and type control
 */
typedef char StlString[80]; /* STL string data requires 80 chars */
typedef unsigned long StlLongInt; /* STL unsigned long integer */
typedef unsigned int StlInt; /* STL unsigned integer */
typedef float StlReal; /* STL requires real data to be float */
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * STL Reader
 */
extern int ReadStlFile(const char *fileName, Polygon *);
/*
 * STL Writer
 */
extern int WriteStlFile(const char *fileName, const Polygon *);
#endif
/* a good practice: end file with a newline */

