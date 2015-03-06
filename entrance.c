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
    while ((1 < argc) && ('-' == argv[1][0])) {
        if (3 > argc) { /* not enough arguments */
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
                if (0 == strcmp(argv[1], "serial")) {
                    control->runMode = 's';
                    break;
                }
                if (0 == strcmp(argv[1], "interact")) {
                    control->runMode = 'i';
                    break;
                }
                if (0 == strcmp(argv[1], "threaded")) {
                    control->runMode = 't';
                    break;
                }
                if (0 == strcmp(argv[1], "mpi")) {
                    control->runMode = 'm';
                    break;
                }
                if (0 == strcmp(argv[1], "gpu")) {
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
    if (1 != argc) {
        fprintf(stderr,"warning, unidentified arguments ignored: %s...\n", argv[1]);
    } 
    /*
     * Check the consistency of run mode and number of processors
     */
    if (('t' == control->runMode) || ('m' == control->runMode) ||
            ('g' == control->runMode)) {
        if (1 >= control->processorN) {
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
    if ('i' == control->runMode) {
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
        if (0 == strncmp(currentLine, "help", sizeof currentLine)) {
            fprintf(stdout, "Options under interactive environment:\n\n");
            fprintf(stdout, "[help]    show this information\n");
            fprintf(stdout, "[init]    generate the initial case input files\n");
            fprintf(stdout, "[solve]   solve current case in serial mode\n");
            fprintf(stdout, "[calc]    access expression calculator\n");
            fprintf(stdout, "[manual]  show a brief user manual of ArtraCFD\n");
            fprintf(stdout, "[exit]    exit program\n");
            continue;
        }
        if (0 == strncmp(currentLine, "init", sizeof currentLine)) {
            GenerateCaseSettingFiles();
            fprintf(stdout, "case files generated successfully\n");
            continue;
        }
        if (0 == strncmp(currentLine, "calc", sizeof currentLine)) {
            ExpressionCalculator();
            continue;
        }
        if (0 == strncmp(currentLine, "manual", sizeof currentLine)) {
            ProgramManual();
            continue;
        }
        if ('\0' == currentLine[0]) { /* no useful information in the command */
            fprintf(stdout, "\n");
            continue;
        }
        if (0 == strncmp(currentLine, "solve", sizeof currentLine)) {
            ShowInformation("Session End");
            return 0;
        }
        if (0 == strncmp(currentLine, "exit", sizeof currentLine)) {
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

