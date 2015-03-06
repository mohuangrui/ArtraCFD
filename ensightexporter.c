/****************************************************************************
 * Functions for Ensight Data Format                                        *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines operations for Ensight gold data format.             *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "ensight.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int WriteEnsightCaseFile(EnsightSet *, const Time *);
static int WriteEnsightGeometryFile(EnsightSet *, const Space *, const Partition *);
static int WriteEnsightVariableFile(const Real *, EnsightSet *, const Space *,
        const Partition *, const Flow *);
static int WriteParticleFile(EnsightSet *, const Particle *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * Ensight transient case file
 * This function initializes an overall transient case file.
 * It's separated because only need to be executed once as initialization.
 */
int InitializeEnsightTransientCaseFile(const Time *time)
{
    ShowInformation("  Initialize Ensight transient case file...");
    FILE *filePointer = fopen("transient.case", "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to ensight case file: transient.case...");
    }
    /* output information to file */
    const EnsightString baseName = "ensight";
    fprintf(filePointer, "FORMAT\n"); 
    fprintf(filePointer, "type: ensight gold\n"); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "GEOMETRY\n"); 
    fprintf(filePointer, "model:            1       %s*****.geo\n", baseName); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "VARIABLE\n"); 
    fprintf(filePointer, "scalar per node:  1  rho  %s*****.rho\n", baseName); 
    fprintf(filePointer, "scalar per node:  1  u    %s*****.u\n", baseName); 
    fprintf(filePointer, "scalar per node:  1  v    %s*****.v\n", baseName); 
    fprintf(filePointer, "scalar per node:  1  w    %s*****.w\n", baseName); 
    fprintf(filePointer, "scalar per node:  1  p    %s*****.p\n", baseName); 
    fprintf(filePointer, "scalar per node:  1  T    %s*****.T\n", baseName); 
    fprintf(filePointer, "vector per node:  1  Vel  %s*****.Vel\n", baseName); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "TIME\n"); 
    fprintf(filePointer, "time set:         1\n"); 
    fprintf(filePointer, "number of steps:          %d\n", (time->totalOutputTimes + 1)); 
    fprintf(filePointer, "filename start number:    0\n"); 
    fprintf(filePointer, "filename increment:       1\n"); 
    fprintf(filePointer, "time values:  "); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * This function write computed data to files with Ensight data format, 
 * including transient and steady output with file names consists of the 
 * default base file name and export step tag. 
 */
int WriteComputedDataEnsight(const Real *U, const Space *space, 
        const Particle *particle, const Time *time, const Partition *part,
        const Flow *flow)
{
    ShowInformation("  Writing field data to file...");
    EnsightSet enSet = { /* initialize Ensight environment */
        .baseName = "ensight", /* data file base name */
        .fileName = {'\0'}, /* data file name */
        .stringData = {'\0'}, /* string data recorder */
    };
    WriteEnsightCaseFile(&enSet, time);
    WriteEnsightGeometryFile(&enSet, space, part);
    WriteEnsightVariableFile(U, &enSet, space, part, flow);
    WriteParticleFile(&enSet, particle);
    return 0;
}
/*
 * Ensight case files
 */
static int WriteEnsightCaseFile(EnsightSet *enSet, const Time *time)
{
    /*
     * Write the steady case file of current step.
     * To get the target file name, the function snprintf is used. It prints 
     * specified format to a string (char array) with buffer overflow 
     * protection. 
     * NOTE: if memeory locations of input objects overlap, the behavior of
     * snprintf is undefined!
     * To make all step number has the same digit number, in the
     * format specifiers, zero padding is added to the integer value. a zero
     * '0' character indicating that zero-padding should be used rather than
     * blank-padding.
     */
    /* store updated basename in filename */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s%05d", enSet->baseName, time->outputCount); 
    /* basename is updated here! */
    snprintf(enSet->baseName, sizeof(EnsightString), "%s", enSet->fileName); 
    /* current filename */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.case", enSet->baseName); 
    FILE *filePointer = fopen(enSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("failed to write data to ensight case file: ensight.case***...");
    }
    /* output information to file */
    fprintf(filePointer, "FORMAT\n"); 
    fprintf(filePointer, "type: ensight gold\n"); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "GEOMETRY\n"); 
    fprintf(filePointer, "model:  %s.geo\n", enSet->baseName); 
    fprintf(filePointer, "\n"); 
    fprintf(filePointer, "VARIABLE\n"); 
    fprintf(filePointer, "constant per case:  Order %d\n", time->outputCount);
    fprintf(filePointer, "constant per case:  Time  %.6g\n", time->currentTime);
    fprintf(filePointer, "constant per case:  Step  %d\n", time->stepCount);
    fprintf(filePointer, "scalar per node:    rho   %s.rho\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    u     %s.u\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    v     %s.v\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    w     %s.w\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    p     %s.p\n", enSet->baseName); 
    fprintf(filePointer, "scalar per node:    T     %s.T\n", enSet->baseName); 
    fprintf(filePointer, "vector per node:    Vel   %s.Vel\n", enSet->baseName); 
    fprintf(filePointer, "\n"); 
    fclose(filePointer); /* close current opened file */
    /*
     * Add the time flag of current export to the transient case
     */
    filePointer = fopen("transient.case", "a");
    if (NULL == filePointer) {
        FatalError("failed to add data to ensight case file: transient.case...");
    }
    if ((time->outputCount % 5) == 0) { /* print to a new line every x outputs */
        fprintf(filePointer, "\n"); 
    }
    fprintf(filePointer, "%.6g ", time->currentTime); 
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * Ensight geometry file
 */
static int WriteEnsightGeometryFile(EnsightSet *enSet, const Space *space, const Partition *part)
{
    /*
     * Write the geometry file (Binary Form).
     * Ensight Maximums: maximum number of nodes in a part is 2GB.
     */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.geo", enSet->baseName);
    FILE *filePointer = fopen(enSet->fileName, "wb");
    if (NULL == filePointer) {
        FatalError("failed to write data to ensight geometry file: ensight.geo***...");
    }
    /*
     * Output information to file, need to strictly follow the Ensight data format.
     * NOTE: if memeory locations of input objects overlap, the behavior of
     * strncpy is undefined!
     * In fwrite, the first size is the sizeof an object, which is given in the
     * units of chars, And the second size (count) is the number of object 
     * that need to be written.
     */
    /* description  at the beginning */
    strncpy(enSet->stringData, "C Binary", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "Ensight Geometry File", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "Written by ArtraCFD", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    /* node id and extents settings */
    strncpy(enSet->stringData, "node id off", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    strncpy(enSet->stringData, "element id off", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    /*
     * Begin to write each part
     */
    int idx = 0; /* linear array index math variable */
    int nodeCount[3] = {0, 0, 0}; /* i j k node number in each part */
    int blankID = 0; /* Ensight geometry iblank entry */
    EnsightReal data = 0; /* the ensight data format */
    for (int partCount = 0, partNum = 1; partCount < part->subN; ++partCount, ++partNum) {
        strncpy(enSet->stringData, "part", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        fwrite(&partNum, sizeof(int), 1, filePointer);
        strncpy(enSet->stringData, part->name[partCount], sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        strncpy(enSet->stringData, "block iblanked", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        /* this line is the total number of nodes in i, j, k */
        nodeCount[0] = (part->iSup[partCount] - part->iSub[partCount]); 
        nodeCount[1] = (part->jSup[partCount] - part->jSub[partCount]);
        nodeCount[2] = (part->kSup[partCount] - part->kSub[partCount]);
        fwrite(nodeCount, sizeof(int), 3, filePointer);
        /* now output the x coordinates of all nodes in current part */
        for (int k = part->kSub[partCount]; k < part->kSup[partCount]; ++k) {
            for (int j = part->jSub[partCount]; j < part->jSup[partCount]; ++j) {
                for (int i = part->iSub[partCount]; i < part->iSup[partCount]; ++i) {
                    data = space->xMin + (i - space->ng) * space->dx;
                    fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                }
            }
        }
        /* now output the y coordinates of all nodes in current part */
        for (int k = part->kSub[partCount]; k < part->kSup[partCount]; ++k) {
            for (int j = part->jSub[partCount]; j < part->jSup[partCount]; ++j) {
                for (int i = part->iSub[partCount]; i < part->iSup[partCount]; ++i) {
                    data = space->yMin + (j - space->ng) * space->dy;
                    fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                }
            }
        }
        /* now output the z coordinates of all nodes in current part */
        for (int k = part->kSub[partCount]; k < part->kSup[partCount]; ++k) {
            for (int j = part->jSub[partCount]; j < part->jSup[partCount]; ++j) {
                for (int i = part->iSub[partCount]; i < part->iSup[partCount]; ++i) {
                    data = space->zMin + (k - space->ng) * space->dz;
                    fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                }
            }
        }
        /*
         * Now output the iblanked array of all nodes in current part
         * blankID = 1 is interior type, will be created and used.
         * blankID = 0 is exterior type, they are blanked-out nodes 
         * and will not be created in the geometry.
         * blankID > 1 or < 0 are any kind of boundary nodes.
         * Transforming from the nodeFlag to blankID are required.
         */
        for (int k = part->kSub[partCount]; k < part->kSup[partCount]; ++k) {
            for (int j = part->jSub[partCount]; j < part->jSup[partCount]; ++j) {
                for (int i = part->iSub[partCount]; i < part->iSup[partCount]; ++i) {
                    idx = (k * space->jMax + j) * space->iMax + i;
                    if (-10 >= space->nodeFlag[idx]) { /* solid region */
                        blankID = 0;
                    } else {
                        if (0 == space->nodeFlag[idx]) { /* fluid */
                            blankID = 1;
                        } else {
                            blankID = space->nodeFlag[idx];
                        }
                    }
                    fwrite(&(blankID), sizeof(int), 1, filePointer);
                }
            }
        }
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * Ensight variables files
 * The values for each node of the structured block are output in 
 * the same IJK order as the coordinates. (The number of nodes in the
 * part are obtained from the corresponding EnSight Gold geometry file.)
 */
static int WriteEnsightVariableFile(const Real *U, EnsightSet *enSet,
        const Space *space, const Partition *part, const Flow *flow)
{
    FILE *filePointer = NULL;
    int idx = 0; /* linear array index math variable */
    EnsightReal data = 0; /* the ensight data format */
    /*
     * Write the scalar field (Binary Form)
     * There are six primitive variables need to be written:
     * density, u, v, w, pressure, temperature
     */
    /*
     * Define the primitive field variables.
     */
    Real rho = 0; 
    Real u = 0;
    Real v = 0;
    Real w = 0;
    Real eT = 0;
    const char nameSuffix[6][5] = {"rho", "u", "v", "w", "p", "T"};
    for (int dim = 0; dim < 6; ++dim) {
        snprintf(enSet->fileName, sizeof(EnsightString), "%s.%s", enSet->baseName, nameSuffix[dim]);
        filePointer = fopen(enSet->fileName, "wb");
        if (NULL == filePointer) {
            FatalError("failed to write data to ensight data file: ensight.***...");
        }
        /* first line description per file */
        strncpy(enSet->stringData, "scalar variable", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        for (int partCount = 0, partNum = 1; partCount < part->subN; ++partCount, ++partNum) {
            /* binary file format */
            strncpy(enSet->stringData, "part", sizeof(EnsightString));
            fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            fwrite(&partNum, sizeof(int), 1, filePointer);
            strncpy(enSet->stringData, "block", sizeof(EnsightString));
            fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
            /* now output the scalar value at each node in current part */
            for (int k = part->kSub[partCount]; k < part->kSup[partCount]; ++k) {
                for (int j = part->jSub[partCount]; j < part->jSup[partCount]; ++j) {
                    for (int i = part->iSub[partCount]; i < part->iSup[partCount]; ++i) {
                        idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                        rho = U[idx];
                        switch (dim) {
                            case 0: /* rho */
                                data = rho;
                                break;
                            case 1: /* u */
                                data = U[idx+1] / rho;
                                break;
                            case 2: /* v */
                                data = U[idx+2] / rho;
                                break;
                            case 3: /* w */
                                data = U[idx+3] / rho;
                                break;
                            case 4: /* p */
                                u = U[idx+1] / rho;
                                v = U[idx+2] / rho;
                                w = U[idx+3] / rho;
                                eT = U[idx+4] / rho;
                                data = (flow->gamma - 1) * rho * (eT - 0.5 * (u * u + v * v + w * w));
                                break;
                            case 5: /* T */
                                u = U[idx+1] / rho;
                                v = U[idx+2] / rho;
                                w = U[idx+3] / rho;
                                eT = U[idx+4] / rho;
                                data = (eT - 0.5 * (u * u + v * v + w * w)) / flow->cv;
                                break;
                            default:
                                break;
                        }
                        fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                    }
                }
            }
        }
        fclose(filePointer); /* close current opened file */
    }
    /*
     * Write the velocity vector field (Binary Form)
     */
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.Vel", enSet->baseName);
    filePointer = fopen(enSet->fileName, "wb");
    if (NULL == filePointer) {
        FatalError("failed to write ensight data file: ensight.Vel***...");
    }
    /* binary file format */
    strncpy(enSet->stringData, "vector variable", sizeof(EnsightString));
    fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
    for (int partCount = 0, partNum = 1; partCount < part->subN; ++partCount, ++partNum) {
        strncpy(enSet->stringData, "part", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        fwrite(&partNum, sizeof(int), 1, filePointer);
        strncpy(enSet->stringData, "block", sizeof(EnsightString));
        fwrite(enSet->stringData, sizeof(char), sizeof(EnsightString), filePointer);
        /*
         * Now output the vector components at each node in current part
         * dimension index of u, v, w is 1, 2, 3 in U in each part, 
         * write u, v, w sequentially
         */
        for (int dim = 1; dim < 4; ++dim) {
            for (int k = part->kSub[partCount]; k < part->kSup[partCount]; ++k) {
                for (int j = part->jSub[partCount]; j < part->jSup[partCount]; ++j) {
                    for (int i = part->iSub[partCount]; i < part->iSup[partCount]; ++i) {
                        idx = ((k * space->jMax + j) * space->iMax + i) * 5;
                        rho = U[idx];
                        switch (dim) {
                            case 1: /* u */
                                data = U[idx+1] / rho;
                                break;
                            case 2: /* v */
                                data = U[idx+2] / rho;
                                break;
                            case 3: /* w */
                                data = U[idx+3] / rho;
                                break;
                            default:
                                break;
                        }
                        fwrite(&data, sizeof(EnsightReal), 1, filePointer);
                    }
                }
            }
        }
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * This file stores the particle information, this will not processed by
 * Ensight, but used for restart.
 */
static int WriteParticleFile(EnsightSet *enSet, const Particle *particle)
{
    snprintf(enSet->fileName, sizeof(EnsightString), "%s.particle", enSet->baseName);
    FILE *filePointer = fopen(enSet->fileName, "w");
    if (NULL == filePointer) {
        FatalError("faild to write particle data file: ensight.particle***...");
    }
    fprintf(filePointer, "N: %d\n", particle->totalN); /* number of objects */
    for (int geoCount = 0; geoCount < particle->totalN; ++geoCount) {
        fprintf(filePointer, "%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g\n", 
                particle->x[geoCount], particle->y[geoCount],
                particle->z[geoCount], particle->r[geoCount],
                particle->u[geoCount], particle->v[geoCount],
                particle->w[geoCount]);
    }
    fclose(filePointer); /* close current opened file */
    return 0;
}
/* a good practice: end file with a newline */

