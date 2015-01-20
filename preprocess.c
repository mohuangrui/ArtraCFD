/****************************************************************************
 * Preprocess                                                               *
 * Last-modified: 20 Jan 2015 12:47:52 AM
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file functions as a preprocessor                                  *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "preprocess.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
#include "casefilegenerator.h"
#include "calculator.h"
#include "casedataloader.h"
#include "casedatachecker.h"
#include "cfdparameters.h"
#include "domainpartition.h"
#include "geometryloader.h"
#include "gcibm.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int Preamble(void);
static int ProgramMemoryAllocate(Field *, Flux *, Space *);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * This is the overall preprocessing function
 */
int Preprocess(Field *field, Flux *flux, Space *space, Particle *particle,
        Time *time, Partition *part, Fluid *fluid, Flow *flow, Reference *reference)
{
    Preamble();
    LoadCaseSettingData(space, time, fluid, reference);
    CheckCaseSettingData(space, time, fluid, reference);
    ComputeCFDParameters(space, time, fluid, flow, reference);
    DomainPartition(part, space);
    ProgramMemoryAllocate(field, flux, space);
    LoadGeometryData(particle, time);
    InitializeDomainGeometryGCIBM(space, particle, part);
    return 0;
}
/*
 * This function prints some general information of the program.
 */
static int Preamble(void)
{
    printf("**********************************************************\n");
    printf("*                        ArtraCFD                        *\n");
    printf("*                     By Huangrui Mo                     *\n");
    printf("**********************************************************\n\n");
    printf("Enter 'help' for a brief user manual\n");
    printf("**********************************************************\n\n");
    char currentLine[200] = {'\0'}; /* store the current read line */
    while (1) {
        printf("\nArtraCFD << ");
        fgets(currentLine, sizeof currentLine, stdin); /* read a line */
        CommandLineProcessor(currentLine); /* process current line */
        printf("\n");
        if (strncmp(currentLine, "help", sizeof currentLine) == 0) {
            printf("Operation options:\n");
            printf("[help]    show this information\n");
            printf("[init]    generate case input files\n");
            printf("[solve]   solve the case\n");
            printf("[calc]    access expression calculator\n");
            printf("[exit]    exit program\n");
            continue;
        }
        if (strncmp(currentLine, "init", sizeof currentLine) == 0) {
            GenerateCaseSettingFiles();
            printf("case files generated successfully\n");
            continue;
        }
        if (strncmp(currentLine, "calc", sizeof currentLine) == 0) {
            ExpressionCalculator();
            continue;
        }
        if (currentLine[0] == '\0') { /* no useful information in the command */
            printf("\n");
            continue;
        }
        if (strncmp(currentLine, "solve", sizeof currentLine) == 0) {
            ShowInformation("Session End");
            return 0;
        }
        if (strncmp(currentLine, "exit", sizeof currentLine) == 0) {
            ShowInformation("Session End");
            exit(EXIT_SUCCESS);
        } 
        /* if non of above is true, then unknow commands */
        printf("artracfd: %s: command not found\n", currentLine);
    }
}
/*
 * This function together with some subfuctions realize the dynamic
 * memory allocation for each global pointer.
 */
static int ProgramMemoryAllocate(Field *field, Flux *flux, Space *space)
{
    ShowInformation("Allocating memory for program...");
    /*
     * Conservative flow variables: rho, rho_u, rho_v, rho_w, E
     */
    int dimU = 5; /* dimension of field variable U */
    int idxMax = dimU * space->kMax * space->jMax * space->iMax;
    field->U = AssignStorage(idxMax, "double");
    field->Un = AssignStorage(idxMax, "double");
    field->Um = AssignStorage(idxMax, "double");
    flux->Fx = AssignStorage(idxMax, "double");
    flux->Fy = AssignStorage(idxMax, "double");
    flux->Fz = AssignStorage(idxMax, "double");
    flux->Gx = AssignStorage(idxMax, "double");
    flux->Gy = AssignStorage(idxMax, "double");
    flux->Gz = AssignStorage(idxMax, "double");
    /*
     * Primitive flow variables: rho, u, v, w, p, T
     */
    idxMax = (dimU + 1) * space->kMax * space->jMax * space->iMax;
    field->Uo = AssignStorage(idxMax, "double");
    /*
     * Node type identifier
     */
    idxMax = space->kMax * space->jMax * space->iMax; /* max index */
    space->ghostFlag = AssignStorage(idxMax, "int");
    space->geoID = AssignStorage(idxMax, "int");
    ShowInformation("Session End");
    return 0;
}
/* a good practice: end file with a newline */

