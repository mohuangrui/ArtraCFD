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
#ifndef ARTRACFD_ENSIGHT_H_ /* if this is the first definition */
#define ARTRACFD_ENSIGHT_H_ /* a unique marker for this header file */
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * Ensight data format and type control
 */
typedef char EnsightString[80]; /* Ensight string data requires 80 chars */
typedef float EnsightReal; /* Ensight requires real data to be float */
/*
 * Ensight configuration structure
 */
typedef struct {
    EnsightString rootName; /* data file root name */
    EnsightString baseName; /* data file base name */
    EnsightString fileName; /* store current open file name */
    EnsightString stringData; /* Ensight string data */
}EnsightSet;
#endif
/* a good practice: end file with a newline */

