/****************************************************************************
 * Geometry Data Loader                                                     *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a loader for geometry data                           *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "geometrystream.h"
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
    if (0 == time->restart) { /* if non-restart, read input file */
        NonrestartGeometryLoader(particle);
    } else { /* if restart, read the particle file */
        RestartGeometryLoader(particle);
    }
    return 0;
}
static int NonrestartGeometryLoader(Particle *particle)
{
    ShowInformation("Loading geometry data ...");
    FILE *filePointer = fopen("artracfd.geo", "r");
    if (NULL == filePointer) {
        FatalError("failed to open geometry file: artracfd.geo...");
    }
    /* read and process file line by line */
    char currentLine[200] = {'\0'}; /* store the current read line */
    int entryCount = 0; /* entry count */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, "sphere count begin", sizeof currentLine)) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(particle->totalN)); 
            continue;
        }
        if (0 == strncmp(currentLine, "sphere begin", sizeof currentLine)) {
            ++entryCount;
            ReadGeometryData(&filePointer, particle);
        }
        continue;
    }
    fclose(filePointer); /* close current opened file */
    /* Check missing information section in configuration */
    if (2 != entryCount) {
        FatalError("missing or repeated necessary information section");
    }
    ShowInformation("Session End");
    return 0;
}
static int RestartGeometryLoader(Particle *particle)
{
    ShowInformation("Restore geometry data ...");
    /* restore geometry information from particle file. */
    FILE *filePointer = fopen("restart.ptk", "r");
    if (NULL == filePointer) {
        FatalError("failed to open restart particle file: restart.ptk...");
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
    if (0 == particle->totalN) { /* no internal geometries */
        return 0;
    }
    /* 
     * Assign storage to store particle information:
     * x, y, z, r, density, u, v, w,       fx, fy, fz, tally 
     * 0, 1, 2, 3,    4,    5, 6, 7,        8,  9, 10,  11      total: 12
     *    need to be read in                  calculated
     */
    particle->entryN = 12;
    particle->headAddress = AssignStorage(particle->totalN * particle->entryN, "Real");
    /* read and store data per object*/
    char currentLine[500] = {'\0'}; /* store the current read line */
    /* set format specifier according to the type of Real */
    char format[100] = "%lg, %lg, %lg, %lg, %lg, %lg, %lg, %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(format, "%g, %g, %g, %g, %g, %g, %g, %g", sizeof format); /* float type */
    }
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        Real *ptk = particle->headAddress + geoCount * particle->entryN;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, format, ptk + 0, ptk + 1, ptk + 2, ptk + 3, ptk + 4, ptk + 5,
                ptk + 6, ptk + 7);
    }
    *filePointerPointer = filePointer; /* return a updated value of file pointer */
    return 0;
}
/*
 * Write geometry information for restart.
 */
int WriteGeometryData(const Particle *particle, const Time *time)
{
    char fileName[25] = {'\0'};
    snprintf(fileName, sizeof(fileName), "restart%05d.ptk", time->outputCount);
    FILE *filePointer = fopen(fileName, "w");
    if (NULL == filePointer) {
        FatalError("faild to write particle data file...");
    }
    fprintf(filePointer, "N: %d\n", particle->totalN); /* number of objects */
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        const Real *ptk = particle->headAddress + geoCount * particle->entryN;
        fprintf(filePointer, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n",
                ptk[0], ptk[1], ptk[2], ptk[3], ptk[4], ptk[5], 
                ptk[6], ptk[7]);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

