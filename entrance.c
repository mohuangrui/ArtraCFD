/****************************************************************************
 * Program Entrance                                                         *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file define functions handle the entering of the program          *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "entrance.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <string.h> /* manipulating strings */
#include "calculator.h"
#include "casefilegenerator.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ConfigureProgram(const Control *);
static int Preamble(void);
static int ProgramManual(void);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
/*
 * This is the overall entrance function
 */
int ProgramEntrance(int argc, char *argv[], Control *control)
{
    /*
     * Loop for command line options
     */
    while ((argc > 1) && (argv[1][0] == '-')) {
        if (argc < 3) { /* not enough arguments */
            fprintf(stderr,"error, empty entry after %s\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        /*
         * argv[1][1] is the actual option character
         */
        switch (argv[1][1]) {
            /*
             * run mode: -m [serial], [interact], [threaded], [mpi], [gpu]
             */
            case 'm':
                ++argv;
                --argc;
                if (strcmp(argv[1], "serial") == 0) {
                    control->runMode = 's';
                    break;
                }
                if (strcmp(argv[1], "interact") == 0) {
                    control->runMode = 'i';
                    break;
                }
                if (strcmp(argv[1], "threaded") == 0) {
                    control->runMode = 't';
                    break;
                }
                if (strcmp(argv[1], "mpi") == 0) {
                    control->runMode = 'm';
                    break;
                }
                if (strcmp(argv[1], "gpu") == 0) {
                    control->runMode = 'g';
                    break;
                }
                fprintf(stderr,"error, bad option %s\n", argv[1]);
                exit(EXIT_FAILURE);
                /*
                 * number of processors: -n N
                 */
            case 'n':
                ++argv;
                --argc;
                sscanf(argv[1], "%d", &(control->processorN));
                break;
            default: 
                fprintf(stderr,"error, bad option %s\n", argv[1]);
                exit(EXIT_FAILURE);
        }
        /*
         * move the argument list up one, move the count down one
         */
        ++argv;
        --argc;
    }
    /*
     * At this point, all the options have been processed.
     * Check to see whether some other information left.
     */
    if (argc != 1) {
        fprintf(stderr,"warning, unidentified arguments ignored: %s...\n", argv[1]);
    } 
    /*
     * Check the consistency of run mode and number of processors
     */
    if ((control->runMode == 't') || (control->runMode == 'm') ||
            (control->runMode == 'g')) {
        if (control->processorN <= 1) {
            fprintf(stderr,"error, illegal number of processors %d\n", control->processorN);
            exit(EXIT_FAILURE);
        }
    } else {
        control->processorN = 1; /* otherwise always set to 1 and ignore input */
    }
    /*
     * Configure program according to inputted commands and information
     */
    ConfigureProgram(control);
    return 0;
}
static int ConfigureProgram(const Control *control)
{
    if (control->runMode == 'i') {
        Preamble();
    }
    return 0;
}
/*
 * This function provides the interactive environment of the program
 */
static int Preamble(void)
{
    fprintf(stdout, "**********************************************************\n");
    fprintf(stdout, "*                        ArtraCFD                        *\n");
    fprintf(stdout, "*                     By Huangrui Mo                     *\n");
    fprintf(stdout, "**********************************************************\n\n");
    fprintf(stdout, "Enter 'help' for more information\n");
    fprintf(stdout, "**********************************************************\n\n");
    char currentLine[200] = {'\0'}; /* store the current read line */
    while (1) {
        fprintf(stdout, "\nArtraCFD << ");
        fgets(currentLine, sizeof currentLine, stdin); /* read a line */
        CommandLineProcessor(currentLine); /* process current line */
        fprintf(stdout, "\n");
        if (strncmp(currentLine, "help", sizeof currentLine) == 0) {
            fprintf(stdout, "Options under interactive environment:\n\n");
            fprintf(stdout, "[help]    show this information\n");
            fprintf(stdout, "[init]    generate case input files\n");
            fprintf(stdout, "[solve]   solve the case in serial mode\n");
            fprintf(stdout, "[calc]    access expression calculator\n");
            fprintf(stdout, "[manual]  show a brief user manual of ArtraCFD\n");
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
        if (strncmp(currentLine, "manual", sizeof currentLine) == 0) {
            ProgramManual();
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
static int ProgramManual(void)
{
    fprintf(stdout, "\n            ArtraCFD User Manual\n\n");
    fprintf(stdout, "SYSNOPSIS:\n");
    fprintf(stdout, "        artracfd [-m runmode] [-n nprocessors]\n");
    fprintf(stdout, "OPTIONS:\n");
    fprintf(stdout, "        -m runmode        run mode: serial, interact, threaded, mpi, gpu\n");
    fprintf(stdout, "        -n nprocessors    number of processors\n");
    fprintf(stdout, "NOTES:\n");
    fprintf(stdout, "        default run mode is 'interact'\n");
    return 0;
}
/* a good practice: end file with a newline */

