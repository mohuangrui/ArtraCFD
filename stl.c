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
 * Required Header Files
 ****************************************************************************/
#include "stl.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include <stdint.h> /* fixed width integer types */
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
typedef enum {
    STLSTR = 80, /* STL header characters */
} StlConst;
/*
 * STL data format and type control
 */
typedef uint8_t StlChar; /* STL char */
typedef StlChar StlStr[STLSTR]; /* STL string */
typedef uint16_t StlInt; /* STL unsigned integer */
typedef uint32_t StlLint; /* STL unsigned long integer */
typedef float StlReal; /* STL real data */
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
/****************************************************************************
 * Function definitions
 ****************************************************************************/
void ReadStlFile(const char *fname, Polyhedron *poly)
{
    StlStr header = {'\0'};
    StlLint facetN = 0;
    StlInt attribute = 0;
    StlReal facetData = 0.0;
    FILE *fp = Fopen(fname, "rb");
    Fread(header, sizeof(StlStr), 1, fp);
    Fread(&facetN, sizeof(StlLint), 1, fp);
    poly->faceN = facetN;
    poly->facet = AssignStorage(poly->faceN * sizeof(*poly->facet));
    for (StlLint n = 0; n < facetN; ++n) {
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].N[X] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].N[Y] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].N[Z] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].v0[X] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].v0[Y] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].v0[Z] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].v1[X] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].v1[Y] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].v1[Z] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].v2[X] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].v2[Y] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, fp);
        poly->facet[n].v2[Z] = facetData;
        Fread(&attribute, sizeof(StlInt), 1, fp);
    }
    fclose(fp);
    return;
}
void WriteStlFile(const char *fname, const Polyhedron *poly)
{
    StlStr header = "binary stl";
    StlLint facetN = 0;
    StlInt attribute = 0;
    StlReal facetData = 0.0;
    FILE *fp = Fopen(fname, "wb");
    fwrite(header, sizeof(StlStr), 1, fp);
    facetN = poly->faceN;
    fwrite(&facetN, sizeof(StlLint), 1, fp);
    for (StlLint n = 0; n < facetN; ++n) {
        facetData = poly->facet[n].N[X];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].N[Y];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].N[Z];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].v0[X];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].v0[Y];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].v0[Z];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].v1[X];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].v1[Y];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].v1[Z];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].v2[X];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].v2[Y];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        facetData = poly->facet[n].v2[Z];
        fwrite(&facetData, sizeof(StlReal), 1, fp);
        fwrite(&attribute, sizeof(StlInt), 1, fp);
    }
    fclose(fp);
    return;
}
/* a good practice: end file with a newline */

