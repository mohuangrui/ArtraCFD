/****************************************************************************
 * Case Data Loader                                                         *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a loader for the case setting file artracfd.case     *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "casedataloader.h"
#include <stdio.h> /* standard library for input and output */
#include <string.h> /* manipulating strings */
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ReadCaseSettingData(Space *, Time *, Flow *, Partition *);
static int ReadBoundaryData(FILE **, Partition *, const int);
static int WriteBoundaryData(FILE **, const Partition *, const int);
static int WriteRegionalInitializerData(FILE **, const int, const Real **);
static int WriteVerifyData(const Space *, const Time *, const Flow *, const Partition *);
static int CheckCaseSettingData(const Space *, const Time *, const Flow *, const Partition *);
/****************************************************************************
 * Function definitions
 ****************************************************************************/
/*
 * This function load the case settings from the case file.
 */
int LoadCaseSettingData(Space *space, Time *time, Flow *flow, Partition *part)
{
    ShowInformation("Loading case setting data ...");
    ReadCaseSettingData(space, time, flow, part);
    WriteVerifyData(space, time, flow, part);
    CheckCaseSettingData(space, time, flow, part);
    ShowInformation("Session End");
    return 0;
}
/*
 * This function read the case settings from the case file.
 * The key is to read and process file line by line. Use "*** begin" 
 * in the case file to identify and control the reading. 
 * The function scanf is notorious for its poor end-of-line handling.
 * Instead, use fgets to read a line of input and sscanf to process it.
 * Note: use a large enough number when using fgets to ensure reading
 * a whole line at a time. fgets will get the entire line including
 * the newline character (\n).
 * NOTE: if memory locations of input objects overlap, the behavior of
 * sscanf is undefined!
 * NOTE: sscanf can correctly handle any space in the target string as
 * well as in the format specifier, therefore, no need to process those
 * lines that will be processed by sscanf.
 * Footnote: In fprintf(), the rvalue type promotions are expected. %f and 
 * %g actually correspond to parameters of type double. Thus in fprintf()
 * there is no difference between %f and %lf, or between %g and %lg. However, 
 * in sscanf() what is passed is a pointer to the variable so no rvalue type 
 * promotions occur or are expected. Thus %f and %lf are quite different in
 * sscanf, but the same in fprintf. Consequently, we need to use %g for 
 * double in fprintf and %lg for double in sscanf. It doesn't matter which
 * you use for fprintf because the fprintf library function treats them as
 * synonymous, but it's crucial to get it right for sscanf. 
 */
static int ReadCaseSettingData(Space *space, Time *time, Flow *flow, Partition *part)
{
    FILE *filePointer = fopen("artracfd.case", "r");
    if (filePointer == NULL) {
        FatalError("failed to open case data file: artracfd.case...");
    }
    /*
     * Read file line by line to get case setting data
     */
    char currentLine[200] = {'\0'}; /* store the current read line */
    int entryCount = 0; /* entry count */
    /* set format specifier according to the type of Real */
    char formatI[5] = "%lg"; /* default is double type */
    char formatIII[15] = "%lg, %lg, %lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
        strncpy(formatIII, "%g, %g, %g", sizeof formatIII); /* float type */
    }
    Real *valueIC = part->valueIC; /* auxiliary pointer for reading regional IC */
    while (fgets(currentLine, sizeof currentLine, filePointer) != NULL) {
        CommandLineProcessor(currentLine); /* process current line */
        if (strncmp(currentLine, "space begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, 
                    &(space->xMin), &(space->yMin), &(space->zMin)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, 
                    &(space->xMax), &(space->yMax), &(space->zMax)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d, %d, %d", 
                    &(space->nx), &(space->ny), &(space->nz)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(space->ng)); 
            continue;
        }
        if (strncmp(currentLine, "time begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->restart)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(time->totalTime)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->totalStep)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(time->numCFL)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, "%d", &(time->totalOutputTimes)); 
            continue;
        }
        if (strncmp(currentLine, "fluid begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refPr)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refMu)); 
            continue;
        }
        if (strncmp(currentLine, "reference begin", sizeof currentLine) == 0) {
            ++entryCount;
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refLength)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refDensity)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refVelocity)); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(flow->refTemperature)); 
            continue;
        }
        if (strncmp(currentLine, "initialization begin", sizeof currentLine) == 0) {
            ++entryCount;
            /*
             * Read initial values rho, u, v, w, p to inner part 0
             */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(part->valueBC[0][0])); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(part->valueBC[0][1])); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(part->valueBC[0][2])); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(part->valueBC[0][3])); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, &(part->valueBC[0][4])); 
            continue;
        }
        if (strncmp(currentLine, "west boundary begin", sizeof currentLine) == 0) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 1);
            continue;
        }
        if (strncmp(currentLine, "east boundary begin", sizeof currentLine) == 0) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 2);
            continue;
        }
        if (strncmp(currentLine, "south boundary begin", sizeof currentLine) == 0) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 3);
            continue;
        }
        if (strncmp(currentLine, "north boundary begin", sizeof currentLine) == 0) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 4);
            continue;
        }
        if (strncmp(currentLine, "front boundary begin", sizeof currentLine) == 0) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 5);
            continue;
        }
        if (strncmp(currentLine, "back boundary begin", sizeof currentLine) == 0) {
            ++entryCount;
            /* Read boundary values for inner part */
            ReadBoundaryData(&filePointer, part, 6);
            continue;
        }
        if (strncmp(currentLine, "plane initialization begin", sizeof currentLine) == 0) {
            /* optional entry do not increase entry count */
            ++part->typeIC[0]; /* regional initializer count and pointer */
            part->typeIC[part->typeIC[0]] = 1; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, valueIC + 0, valueIC + 1, valueIC + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, valueIC + 3, valueIC + 4, valueIC + 5); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 6); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 7); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 8); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 9); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 10); 
            valueIC = valueIC + 11; /* update data loader pointer */
            continue;
        }
        if (strncmp(currentLine, "sphere initialization begin", sizeof currentLine) == 0) {
            /* optional entry do not increase entry count */
            ++part->typeIC[0]; /* regional initializer count and pointer */
            part->typeIC[part->typeIC[0]] = 2; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, valueIC + 0, valueIC + 1, valueIC + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 3); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 4); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 5); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 6); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 7); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 8); 
            valueIC = valueIC + 9; /* update data loader pointer */
            continue;
        }
        if (strncmp(currentLine, "box initialization begin", sizeof currentLine) == 0) {
            /* optional entry do not increase entry count */
            ++part->typeIC[0]; /* regional initializer count and pointer */
            part->typeIC[part->typeIC[0]] = 3; /* IC type id */
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, valueIC + 0, valueIC + 1, valueIC + 2); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatIII, valueIC + 3, valueIC + 4, valueIC + 5); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 6); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 7); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 8); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 9); 
            fgets(currentLine, sizeof currentLine, filePointer);
            sscanf(currentLine, formatI, valueIC + 10); 
            valueIC = valueIC + 11; /* update data loader pointer */
            continue;
        }
    }
    fclose(filePointer); /* close current opened file */
    /*
     * Check missing information section in configuration
     */
    if (entryCount != 11) {
        FatalError("missing or repeated necessary information section");
    }
    return 0;
}
/*
 * Boundary type ID:
 * 0: interior region (default value)
 * 1: inlet
 * 2: outflow
 * 3: slip wall
 * 4: nonslip wall
 * 5: periodic
 * -5: periodic pair
 */
static int ReadBoundaryData(FILE **filePointerPointer, Partition *part, const int partID)
{
    if (part->typeBC[partID] == -5) { /* already set as periodic pair */
        return 0;
    }
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    char currentLine[200] = {'\0'}; /* store the current read line */
    char formatI[5] = "%lg"; /* default is double type */
    if (sizeof(Real) == sizeof(float)) { /* if set Real as float */
        strncpy(formatI, "%g", sizeof formatI); /* float type */
    }
    fgets(currentLine, sizeof currentLine, filePointer);
    CommandLineProcessor(currentLine); /* process current line */
    if (strncmp(currentLine, "inlet", sizeof currentLine) == 0) {
        part->typeBC[partID] = 1;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[partID][0])); 
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[partID][1])); 
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[partID][2])); 
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[partID][3])); 
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[partID][4])); 
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    if (strncmp(currentLine, "outflow", sizeof currentLine) == 0) {
        part->typeBC[partID] = 2;
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    if (strncmp(currentLine, "slip wall", sizeof currentLine) == 0) {
        part->typeBC[partID] = 3;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[partID][5])); 
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    if (strncmp(currentLine, "nonslip wall", sizeof currentLine) == 0) {
        part->typeBC[partID] = 4;
        fgets(currentLine, sizeof currentLine, filePointer);
        sscanf(currentLine, formatI, &(part->valueBC[partID][5])); 
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    if (strncmp(currentLine, "periodic", sizeof currentLine) == 0) {
        /* only need to set id and its periodic pair */
        part->typeBC[partID] = 5;
        part->typeBC[partID - part->normalZ[partID] - part->normalY[partID] - part->normalX[partID]] = -5;
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    FatalError("Unidentified boundary type...");
    return 0;
}
static int WriteBoundaryData(FILE **filePointerPointer, const Partition *part, const int partID)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    if (part->typeBC[partID] == 1) {
        fprintf(filePointer, "boundary type: inlet\n"); 
        fprintf(filePointer, "density: %.6g\n", part->valueBC[partID][0]);
        fprintf(filePointer, "x velocity: %.6g\n", part->valueBC[partID][1]);
        fprintf(filePointer, "y velocity: %.6g\n", part->valueBC[partID][2]);
        fprintf(filePointer, "z velocity: %.6g\n", part->valueBC[partID][3]);
        fprintf(filePointer, "pressure: %.6g\n", part->valueBC[partID][4]);
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    if (part->typeBC[partID] == 2) {
        fprintf(filePointer, "boundary type: outflow\n"); 
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    if (part->typeBC[partID] == 3) {
        fprintf(filePointer, "boundary type: slip wall\n"); 
        fprintf(filePointer, "temperature: %.6g\n", part->valueBC[partID][5]);
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    if (part->typeBC[partID] == 4) {
        fprintf(filePointer, "boundary type: nonslip wall\n"); 
        fprintf(filePointer, "temperature: %.6g\n", part->valueBC[partID][5]);
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    if ((part->typeBC[partID] == 5) || (part->typeBC[partID] == -5)) {
        fprintf(filePointer, "boundary type: periodic\n"); 
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        return 0;
    }
    FatalError("Unidentified boundary type...");
    return 0;
}
static int WriteRegionalInitializerData(FILE **filePointerPointer, const int typeIC, const Real **valueICPointerPointer)
{
    FILE *filePointer = *filePointerPointer; /* get the value of file pointer */
    const Real *valueIC = *valueICPointerPointer; /* get the current valueIC pointer */
    if (typeIC == 1) {
        fprintf(filePointer, "regional initialization: plane\n"); 
        fprintf(filePointer, "plane point x, y, z: %.6g, %.6g, %.6g\n", valueIC[0], valueIC[1], valueIC[2]);
        fprintf(filePointer, "plane normal nx, ny, nz: %.6g, %.6g, %.6g\n", valueIC[3], valueIC[4], valueIC[5]);
        fprintf(filePointer, "density: %.6g\n", valueIC[6]);
        fprintf(filePointer, "x velocity: %.6g\n", valueIC[7]);
        fprintf(filePointer, "y velocity: %.6g\n", valueIC[8]);
        fprintf(filePointer, "z velocity: %.6g\n", valueIC[9]);
        fprintf(filePointer, "pressure: %.6g\n", valueIC[10]);
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        *valueICPointerPointer = valueIC + 11; /* update pointer of valueIC queue */
        return 0;
    }
    if (typeIC == 2) {
        fprintf(filePointer, "regional initialization: sphere\n"); 
        fprintf(filePointer, "center point x, y, z: %.6g, %.6g, %.6g\n", valueIC[0], valueIC[1], valueIC[2]);
        fprintf(filePointer, "radius: %.6g\n", valueIC[3]);
        fprintf(filePointer, "density: %.6g\n", valueIC[4]);
        fprintf(filePointer, "x velocity: %.6g\n", valueIC[5]);
        fprintf(filePointer, "y velocity: %.6g\n", valueIC[6]);
        fprintf(filePointer, "z velocity: %.6g\n", valueIC[7]);
        fprintf(filePointer, "pressure: %.6g\n", valueIC[8]);
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        *valueICPointerPointer = valueIC + 9; /* update pointer of valueIC queue */
        return 0;
    }
    if (typeIC == 3) {
        fprintf(filePointer, "regional initialization: box\n"); 
        fprintf(filePointer, "xmin, ymin, zmin: %.6g, %.6g, %.6g\n", valueIC[0], valueIC[1], valueIC[2]);
        fprintf(filePointer, "xmax, ymax, zmax: %.6g, %.6g, %.6g\n", valueIC[3], valueIC[4], valueIC[5]);
        fprintf(filePointer, "density: %.6g\n", valueIC[6]);
        fprintf(filePointer, "x velocity: %.6g\n", valueIC[7]);
        fprintf(filePointer, "y velocity: %.6g\n", valueIC[8]);
        fprintf(filePointer, "z velocity: %.6g\n", valueIC[9]);
        fprintf(filePointer, "pressure: %.6g\n", valueIC[10]);
        *filePointerPointer = filePointer; /* return a updated value of file pointer */
        *valueICPointerPointer = valueIC + 11; /* update pointer of valueIC queue */
        return 0;
    }
    FatalError("Unidentified regional initializer type...");
    return 0;
}
/*
 * This function outputs the case setting data to a file for verification.
 */
static int WriteVerifyData(const Space *space, const Time *time, const Flow *flow, const Partition *part)
{
    ShowInformation("  Data outputted into artracfd.verify...");
    FILE *filePointer = fopen("artracfd.verify", "w");
    if (filePointer == NULL) {
        FatalError("failed to write data to file: artracfd.verify");
    }
    /* output information to file */
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#                     Case Conformation for ArtraCFD                          -\n");
    fprintf(filePointer, "#                                                                             -\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Space Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "domain xmin, ymin, zmin: %.6g, %.6g, %.6g\n", space->xMin, space->yMin, space->zMin); 
    fprintf(filePointer, "domain xmax, ymax, zmax: %.6g, %.6g, %.6g\n", space->xMax, space->yMax, space->zMax); 
    fprintf(filePointer, "x, y, z mesh number: %d, %d, %d\n", space->nx, space->ny, space->nz); 
    fprintf(filePointer, "exterior ghost cell layers: %d\n", space->ng); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                          >> Time Domain <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "restart flag: %d\n", time->restart); 
    fprintf(filePointer, "total evolution time: %.6g\n", time->totalTime); 
    fprintf(filePointer, "maximum number of steps: %d\n", time->totalStep); 
    fprintf(filePointer, "CFL condition number: %.6g\n", time->numCFL); 
    fprintf(filePointer, "exporting data times: %d\n", time->totalOutputTimes); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                    >> Fluid and Flow Properties <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "Prandtl number: %.6g\n", flow->refPr); 
    fprintf(filePointer, "modify coefficient of dynamic viscosity: %.6g\n", flow->refMu); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                        >> Reference Values  <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "length: %.6g\n", flow->refLength); 
    fprintf(filePointer, "density: %.6g\n", flow->refDensity); 
    fprintf(filePointer, "velocity: %.6g\n", flow->refVelocity); 
    fprintf(filePointer, "temperature: %.6g\n", flow->refTemperature); 
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                            >> NOTE <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "# Values in following parts are relative to reference values. Hence, they need\n");
    fprintf(filePointer, "# to be normalized by the given reference values. Like pressure should be\n");
    fprintf(filePointer, "# normalized by reference density times reference velocity square.\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                     >> Flow Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "density: %.6g\n", part->valueBC[0][0]);
    fprintf(filePointer, "x velocity: %.6g\n", part->valueBC[0][1]);
    fprintf(filePointer, "y velocity: %.6g\n", part->valueBC[0][2]);
    fprintf(filePointer, "z velocity: %.6g\n", part->valueBC[0][3]);
    fprintf(filePointer, "pressure: %.6g\n", part->valueBC[0][4]);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                     >> Boundary Condition <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian West\n"); 
    WriteBoundaryData(&filePointer, part, 1);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian East\n"); 
    WriteBoundaryData(&filePointer, part, 2);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian South\n"); 
    WriteBoundaryData(&filePointer, part, 3);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian North\n"); 
    WriteBoundaryData(&filePointer, part, 4);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Front\n"); 
    WriteBoundaryData(&filePointer, part, 5);
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "Domian Back\n"); 
    WriteBoundaryData(&filePointer, part, 6);
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#                  >> Regional Initialization <<\n");
    fprintf(filePointer, "#\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    const Real *valueIC = part->valueIC; /* pointer to the value queue of initial values */
    int i = 0;
    for (i = 1; i <= part->typeIC[0]; ++i) {
        fprintf(filePointer, "#\n");
        WriteRegionalInitializerData(&filePointer, part->typeIC[i], &valueIC);
    }
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fprintf(filePointer, "#------------------------------------------------------------------------------\n");
    fclose(filePointer); /* close current opened file */
    return 0;
}
/*
 * This function do some parameter checking
 */
static int CheckCaseSettingData(const Space *space, const Time *time, const Flow *flow, const Partition *part)
{
    ShowInformation("  Preliminary case data checking ...");
    /* space */
    if (((space->xMax - space->xMin) < 0) || ((space->yMax - space->yMin) < 0) ||
            ((space->zMax - space->zMin) < 0)) {
        FatalError("wrong domian region values in case settings");
    }
    if ((space->nz < 1) || (space->ny < 1) || (space->nx < 1)
            || (space->ng < 1)) {
        FatalError("too small mesh values in case settings");
    }
    /* time */
    if ((time->restart < 0) || (time->restart > 1)|| (time->totalTime <= 0)
            || (time->numCFL <= 0) || (time->totalOutputTimes < 1)
            || (time->totalStep < 1)) {
        FatalError("wrong values in time section of case settings");
    }
    /* fluid and flow */
    if ((flow->refPr <= 0) || (flow->refMu < 0)) {
        FatalError("wrong values in fluid and flow section of case settings");
    }
    /* reference */
    if ((flow->refLength <= 0) || (flow->refDensity <= 0) || 
            (flow->refVelocity <= 0) || (flow->refTemperature <= 0)) {
        FatalError("wrong values in reference section of case settings");
    }
    /* initialization */
    if ((part->valueBC[0][0] < 0) || (part->valueBC[0][1] < 0) || 
            (part->valueBC[0][2] < 0) || (part->valueBC[0][3] < 0) ||
            (part->valueBC[0][4] < 0)) {
        FatalError("wrong values in initialization section of case settings");
    }
    return 0;
}
/* a good practice: end file with a newline */

