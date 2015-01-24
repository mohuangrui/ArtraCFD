/****************************************************************************
 * Preprocess                                                               *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines a functions as a preprocessor                        *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "preprocess.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <string.h> /* manipulating strings */
#include "calculator.h"
#include "casefilegenerator.h"
#include "casedataloader.h"
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
    fprintf(stdout, "**********************************************************\n");
    fprintf(stdout, "*                        ArtraCFD                        *\n");
    fprintf(stdout, "*                     By Huangrui Mo                     *\n");
    fprintf(stdout, "**********************************************************\n\n");
    fprintf(stdout, "Enter 'help' for a brief user manual\n");
    fprintf(stdout, "**********************************************************\n\n");
    char currentLine[200] = {'\0'}; /* store the current read line */
    while (1) {
        fprintf(stdout, "\nArtraCFD << ");
        fgets(currentLine, sizeof currentLine, stdin); /* read a line */
        CommandLineProcessor(currentLine); /* process current line */
        fprintf(stdout, "\n");
        if (strncmp(currentLine, "help", sizeof currentLine) == 0) {
            fprintf(stdout, "Operation options:\n");
            fprintf(stdout, "[help]    show this information\n");
            fprintf(stdout, "[init]    generate case input files\n");
            fprintf(stdout, "[solve]   solve the case\n");
            fprintf(stdout, "[calc]    access expression calculator\n");
            fprintf(stdout, "[exit]    exit program\n");
            continue;
        }
        if (strncmp(currentLine, "init", sizeof currentLine) == 0) {
            GenerateCaseSettingFiles();
            fprintf(stdout, "case files generated successfully\n");
            continue;
        }
        if (strncmp(currentLine, "calc", sizeof currentLine) == 0) {
            ExpressionCalculator();
            continue;
        }
        if (currentLine[0] == '\0') { /* no useful information in the command */
            fprintf(stdout, "\n");
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
        fprintf(stdout, "artracfd: %s: command not found\n", currentLine);
    }
}
/*
 * This function together with some subfuctions realize the dynamic
 * memory allocation for each global pointer. The storage retrieving
 * need to be done in the postprocessor.
 */
static int ProgramMemoryAllocate(Field *field, Flux *flux, Space *space)
{
    ShowInformation("Allocating memory for program...");
    /*
     * Conservative flow variables: rho, rho_u, rho_v, rho_w, E
     */
    int dimU = 5; /* dimension of field variable U */
    int idxMax = dimU * space->kMax * space->jMax * space->iMax;
    field->U = AssignStorage(idxMax, "Real");
    field->Un = AssignStorage(idxMax, "Real");
    field->Um = AssignStorage(idxMax, "Real");
    flux->Fx = AssignStorage(idxMax, "Real");
    flux->Fy = AssignStorage(idxMax, "Real");
    flux->Fz = AssignStorage(idxMax, "Real");
    flux->Gx = AssignStorage(idxMax, "Real");
    flux->Gy = AssignStorage(idxMax, "Real");
    flux->Gz = AssignStorage(idxMax, "Real");
    /*
     * Primitive flow variables: rho, u, v, w, p, T
     */
    idxMax = (dimU + 1) * space->kMax * space->jMax * space->iMax;
    field->Uo = AssignStorage(idxMax, "Real");
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

