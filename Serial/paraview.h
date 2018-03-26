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
#ifndef ARTRACFD_PARAVIEW_H_ /* if this is the first definition */
#define ARTRACFD_PARAVIEW_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * Paraview data format and type control
 */
typedef char ParaviewString[80]; /* string data */
typedef double ParaviewReal; /* real data */
/*
 * Paraview configuration structure
 */
typedef struct {
    ParaviewString rootName; /* data file root name */
    ParaviewString baseName; /* data file base name */
    ParaviewString fileName; /* store current open file name */
    ParaviewString fileExt; /* data file extension */
    ParaviewString intType; /* Paraview int type */
    ParaviewString floatType; /* Paraview float type */
    ParaviewString byteOrder; /* byte order of data */
}ParaviewSet;
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Structured data writer and reader
 */
extern int WriteStructuredDataParaview(const Time *, const Space *, const Model *);
extern int ReadStructuredDataParaview(Time *, Space *, const Model *);
/*
 * Poly data writer and reader
 */
extern int WritePolyDataParaview(const Time *, const Geometry *);
extern int ReadPolyDataParaview(const Time *, Geometry *);
/*
 * Polyhedron status writer and reader
 */
extern int WritePolyhedronStateData(const int start, const int end, FILE *filePointer, const Geometry *);
extern int ReadPolyhedronStateData(const int start, const int end, FILE *filePointer, Geometry *);
#endif
/* a good practice: end file with a newline */

