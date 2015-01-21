/****************************************************************************
 * Geometry Data Loader                                                     *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a loader for geometry data                           *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "geometryloader.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NonrestartGeometryLoader(Particle *);
static int RestartGeometryLoader(Particle *);
static int ReadGeometryData(FILE *, Particle *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function load geometry data from the geometry file.
 */
int LoadGeometryData(Particle *particle, const Time *time)
{
    if (time->restart == 0) { /* if non-restart, read input file */
        NonrestartGeometryLoader(particle);
    } else { /* if restart, read the particle file */
        RestartGeometryLoader(particle);
    }
    return 0;
}
static int NonrestartGeometryLoader(Particle *particle)
{
    ShowInformation("Loading inputed geometry data ...");
    FILE *filePointer = fopen("artracfd.geo", "r");
    if (filePointer == NULL) {
        FatalError("failed to open geometry file");
    }
    /* read and process file line by line */
    char currentLine[200] = {'\0'}; /* store the current read line */
    int entryCount = 0; /* entry count */
    while (fgets(currentLine, sizeof currentLine, filePointer) != NULL) {
        CommandLineProcessor(currentLine); /* process current line */
        if (strncmp(currentLine, "count begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(particle->totalN)); 
            continue;
        }
        if (strncmp(currentLine, "circle begin", sizeof currentLine) == 0) {
            ++entryCount;
            ReadGeometryData(filePointer, particle);
        }
        continue;
    }
    fclose(filePointer); /* close current opened file */
    /* Check missing information section in configuration */
    if (entryCount != 2) {
        FatalError("missing or repeated necessary information section");
    }
    ShowInformation("Session End");
    return 0;
}
static int RestartGeometryLoader(Particle *particle)
{
    ShowInformation("Restore geometry data ...");
    /* restore geometry information from particle file. */
    FILE *filePointer = fopen("restart.particle", "r");
    if (filePointer == NULL) {
        FatalError("failed to open restart particle file");
    }
    /* read and process file line by line */
    char currentLine[200] = {'\0'}; /* store the current read line */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "N: %d", &(particle->totalN)); 
    ReadGeometryData(filePointer, particle);
    fclose(filePointer); /* close current opened file */
    ShowInformation("Session End");
    return 0;
}
static int ReadGeometryData(FILE *filePointer, Particle *particle)
{
    if (particle->totalN == 0) { /* no internal geometries */
        return 0;
    }
    /* first assign storage to particle pointers */
    particle->x = AssignStorage(particle->totalN, "double");
    particle->y = AssignStorage(particle->totalN, "double");
    particle->z = AssignStorage(particle->totalN, "double");
    particle->r = AssignStorage(particle->totalN, "double");
    particle->u = AssignStorage(particle->totalN, "double");
    particle->v = AssignStorage(particle->totalN, "double");
    particle->w = AssignStorage(particle->totalN, "double");
    /* then read and store data per object*/
    char currentLine[200] = {'\0'}; /* store the current read line */
    int geoCount = 0; /* geometry object count */
    for (geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, "%lg, %lg, %lg, %lg, %lg, %lg, %lg", 
                &(particle->x[geoCount]), &(particle->y[geoCount]),
                &(particle->z[geoCount]), &(particle->r[geoCount]),
                &(particle->u[geoCount]), &(particle->v[geoCount]),
                &(particle->w[geoCount]));
    }
    return 0;
}
/* a good practice: end file with a newline */

