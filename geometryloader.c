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
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int NonrestartGeometryLoader(Particle *);
static int RestartGeometryLoader(Particle *);
static int ReadGeometryData(FILE **, Particle *);
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
        FatalError("failed to open geometry file: artracfd.geo...");
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
            ReadGeometryData(&filePointer, particle);
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
        FatalError("failed to open restart particle file: restart.particle...");
    }
    /* read and process file line by line */
    char currentLine[200] = {'\0'}; /* store the current read line */
    fgets(currentLine, sizeof currentLine, filePointer);
    sscanf(currentLine, "N: %d", &(particle->totalN)); 
    ReadGeometryData(&filePointer, particle);
    fclose(filePointer); /* close current opened file */
    ShowInformation("Session End");
    return 0;
}
static int ReadGeometryData(FILE **filePointerPointer, Particle *particle)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    if (particle->totalN == 0) { /* no internal geometries */
        return 0;
    }
    /* first assign storage to particle pointers */
    particle->headAddress = AssignStorage(particle->totalN * 7, "Real");
    particle->x = particle->headAddress + particle->totalN * 0;
    particle->y = particle->headAddress + particle->totalN * 1;
    particle->z = particle->headAddress + particle->totalN * 2;
    particle->r = particle->headAddress + particle->totalN * 3;
    particle->u = particle->headAddress + particle->totalN * 4;
    particle->v = particle->headAddress + particle->totalN * 5;
    particle->w = particle->headAddress + particle->totalN * 6;
    /* then read and store data per object*/
    char currentLine[200] = {'\0'}; /* store the current read line */
    int geoCount = 0; /* geometry object count */
    /* set format specifier according to the type of Real */
    char formatVII[40] = "%lg, %lg, %lg, %lg, %lg, %lg, %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatVII, "%g, %g, %g, %g, %g, %g, %g", sizeof formatVII); /* float type */
    }
    for (geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatVII, 
                &(particle->x[geoCount]), &(particle->y[geoCount]),
                &(particle->z[geoCount]), &(particle->r[geoCount]),
                &(particle->u[geoCount]), &(particle->v[geoCount]),
                &(particle->w[geoCount]));
    }
    *filePointerPointer = filePointer; /* return a updated value of file pointer */
    return 0;
}
/* a good practice: end file with a newline */

