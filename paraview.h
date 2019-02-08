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
#ifndef ARTRACFD_PARAVIEW_H_ /* if undefined */
#define ARTRACFD_PARAVIEW_H_ /* set a unique marker */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    PVSTR = 80, /* string data length */
    PVVARSTR = 10, /* variable name length */
    PVSCAN = 10, /* maximum number of scalar variables */
    PVVECN = 1, /* maximum number of vector variables */
} PvConst;
typedef char PvStr[PVSTR]; /* string data */
typedef Real PvReal; /* real data */
typedef struct {
    PvStr rname; /* data file root name */
    PvStr bname; /* data file base name */
    PvStr fname; /* store current open file name */
    PvStr fext; /* data file extension */
    PvStr fmt; /* format specifier */
    PvStr intType; /* int type */
    PvStr floatType; /* float type */
    PvStr byteOrder; /* byte order of data */
    int scaN; /* number of scalar variables */
    char sca[PVSCAN][PVVARSTR]; /* scalar variables */
    int vecN; /* number of vector variables */
    char vec[PVVECN][PVVARSTR]; /* vector variables */
} PvSet; /* configuration structure */
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Structured data writer and reader
 */
extern void WriteStructuredDataParaview(const Time *, const Space *, const Model *);
extern void ReadStructuredDataParaview(Time *, Space *, const Model *);
/*
 * Poly data writer and reader
 */
extern void WritePolyDataParaview(const Time *, const Geometry *const);
extern void ReadPolyDataParaview(const Time *, Geometry *const);
#endif
/* a good practice: end file with a newline */

