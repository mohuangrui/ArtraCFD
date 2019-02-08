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
#ifndef ARTRACFD_ENSIGHT_H_ /* if undefined */
#define ARTRACFD_ENSIGHT_H_ /* set a unique marker */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    ENSTR = 80, /* string data length */
    ENVARSTR = 10, /* variable name length */
    ENSCAN = 10, /* maximum number of scalar variables */
    ENVECN = 1, /* maximum number of vector variables */
} EnConst;
typedef char EnStr[ENSTR]; /* string data */
typedef float EnReal; /* real data */
typedef struct {
    EnStr rname; /* data file root name */
    EnStr bname; /* data file base name */
    EnStr fname; /* store current open file name */
    EnStr str; /* string data */
    EnStr fmt; /* format specifier */
    EnStr gtag; /* geometry name tag */
    EnStr vtag; /* variable name tag */
    EnStr dtype; /* data type */
    int part[LIMIT]; /* part control */
    int scaN; /* number of scalar variables */
    char sca[ENSCAN][ENVARSTR]; /* scalar variables */
    int vecN; /* number of vector variables */
    char vec[ENVECN][ENVARSTR]; /* vctor variables */
} EnSet; /* configuration structure */
/****************************************************************************
 * Public Functions Declaration
 ****************************************************************************/
/*
 * Structured data writer and reader
 */
extern void WriteStructuredDataEnsight(const Time *, const Space *, const Model *);
extern void ReadStructuredDataEnsight(Time *, Space *, const Model *);
/*
 * Poly data writer and reader
 */
extern void WritePolyDataEnsight(const Time *, const Geometry *const);
extern void ReadPolyDataEnsight(const Time *, Geometry *const);
#endif
/* a good practice: end file with a newline */

