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
#include "commons.h"
/****************************************************************************
 * Data Structure Declarations
 ****************************************************************************/
/*
 * STL data format and type control
 */
typedef char StlString[80]; /* UINT8[80], STL string data requires 80 chars */
typedef unsigned int StlLongInt; /* UINT32, STL unsigned long integer */
typedef unsigned short int StlInt; /* UINT16, STL unsigned integer */
typedef float StlReal; /* REAL32, STL real data */
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ReadStlFile(const char *fileName, Polyhedron *poly)
{
    StlString header = {'\0'};
    StlLongInt facetN = 0;
    StlInt attributeCount = 0;
    StlReal facetData = 0.0;
    FILE *filePointer = fopen(fileName, "rb");
    if (NULL == filePointer) {
        FatalError("failed to open STL file...");
    }
    Fread(header, sizeof(char), sizeof(StlString), filePointer);
    Fread(&facetN, sizeof(StlLongInt), 1, filePointer);
    poly->faceN = facetN;
    poly->facet = AssignStorage(poly->faceN * sizeof(*poly->facet));
    for (StlLongInt n = 0; n < facetN; ++n) {
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].N[X] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].N[Y] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].N[Z] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].v0[X] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].v0[Y] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].v0[Z] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].v1[X] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].v1[Y] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].v1[Z] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].v2[X] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].v2[Y] = facetData;
        Fread(&facetData, sizeof(StlReal), 1, filePointer);
        poly->facet[n].v2[Z] = facetData;
        Fread(&attributeCount, sizeof(StlInt), 1, filePointer);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
int WriteStlFile(const char *fileName, const Polyhedron *poly)
{
    StlString header = {'\0'};
    StlLongInt facetN = 0;
    StlInt attributeCount = 0;
    StlReal facetData = 0.0;
    FILE *filePointer = fopen(fileName, "wb");
    if (NULL == filePointer) {
        FatalError("failed to open STL file...");
    }
    strncpy(header, "binary stl", sizeof(StlString));
    fwrite(header, sizeof(char), sizeof(StlString), filePointer);
    facetN = poly->faceN;
    fwrite(&facetN, sizeof(StlLongInt), 1, filePointer);
    for (StlLongInt n = 0; n < facetN; ++n) {
        facetData = poly->facet[n].N[X];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].N[Y];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].N[Z];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].v0[X];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].v0[Y];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].v0[Z];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].v1[X];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].v1[Y];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].v1[Z];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].v2[X];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].v2[Y];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = poly->facet[n].v2[Z];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        fwrite(&attributeCount, sizeof(StlInt), 1, filePointer);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

