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
 * Required Header Files
 ****************************************************************************/
#include "stl.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
/****************************************************************************
 * Function definitions
 ****************************************************************************/
int ReadStlFile(const char *fileName, Polygon *polygon)
{
    int idx = 0; /* linear array index math variable */
    StlString header = {'\0'};
    StlLongInt facetN = 0;
    StlInt attributeCount = 0;
    StlReal facetData = 0.0;
    FILE *filePointer = fopen(fileName, "rb");
    if (NULL == filePointer) {
        FatalError("failed to read STL file...");
    }
    fread(header, sizeof(char), sizeof(StlString), filePointer);
    fread(&facetN, sizeof(StlLongInt), 1, filePointer);
    polygon->facetN = facetN;
    for (int count = 0; count < facetN; ++count) {
        idx = count * ENTRYFACET;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FNX] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FNY] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FNZ] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FX1] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FY1] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FZ1] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FX2] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FY2] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FZ2] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FX3] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FY3] = facetData;
        fread(&facetData, sizeof(StlReal), 1, filePointer);
        polygon->facet[idx+FZ3] = facetData;
        fread(&attributeCount, sizeof(StlInt), 1, filePointer);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
int WriteStlFile(const char *fileName, const Polygon *polygon)
{
    int idx = 0; /* linear array index math variable */
    StlString header = {'\0'};
    StlLongInt facetN = 0;
    StlInt attributeCount = 0;
    StlReal facetData = 0.0;
    FILE *filePointer = fopen(fileName, "wb");
    if (NULL == filePointer) {
        FatalError("failed to write STL file...");
    }
    strncpy(header, "binary stl", sizeof(StlString));
    fwrite(header, sizeof(char), sizeof(StlString), filePointer);
    facetN = polygon->facetN;
    fwrite(&facetN, sizeof(StlLongInt), 1, filePointer);
    for (int count = 0; count < facetN; ++count) {
        idx = count * ENTRYFACET;
        facetData = polygon->facet[idx+FNX];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FNY];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FNZ];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FX1];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FY1];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FZ1];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FX2];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FY2];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FZ2];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FX3];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FY3];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        facetData = polygon->facet[idx+FZ3];
        fwrite(&facetData, sizeof(StlReal), 1, filePointer);
        fwrite(&attributeCount, sizeof(StlInt), 1, filePointer);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

