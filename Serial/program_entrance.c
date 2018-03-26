/****************************************************************************
 *                              ArtraCFD                                    *
 *                          <By Huangrui Mo>                                *
 * Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>                        *
 * This file is part of ArtraCFD.                                           *
 * ArtraCFD is free software: you can redistribute it and/or modify it      *
 * under the terms of the GNU General Public License as published by        *
 * the Free Software Foundation, either version 3 of the License, or        *
 * (at your option) any later version.                                      *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "program_entrance.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <string.h> /* manipulating strings */
#include "calculator.h"
#include "case_generator.h"
#include "commons.h"
/****************************************************************************
 * Static Function Declarations
 ****************************************************************************/
static int ConfigureProgram(Control *);
static int Preamble(Control *);
static int ProgramManual(void);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
int ProgramEntrance(int argc, char *argv[], Control *control)
{
    int nscan = 0; /* read conversion count */
    /*
     * Loop command line to process options
     */
    while ((1 < argc) && ('-' == argv[1][0])) {
        if (3 > argc) { /* not enough arguments */
            fprintf(stderr,"error, empty entry after %s\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        switch (argv[1][1]) { /* argv[1][1] is the actual option character */
            /*
             * run mode: -m [interact], [serial], [threaded], [mpi], [gpu]
             */
            case 'm':
                ++argv;
                --argc;
                if (0 == strcmp(argv[1], "interact")) {
                    control->runMode = 'i';
                    break;
                }
                if (0 == strcmp(argv[1], "serial")) {
                    control->runMode = 's';
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
                nscan = sscanf(argv[1], "%d", &(control->procN));
                VerifyReadConversion(nscan, 1);
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
     * Configure program according to inputted commands and information
     */
    ConfigureProgram(control);
    return 0;
}
static int ConfigureProgram(Control *control)
{
    switch (control->runMode) {
        case 'i': /* interaction mode */
            Preamble(control);
            break;
        case 's': /* serial mode */
            break;
        case 't': /* threaded mode */
            break;
        case 'm': /* mpi mode */
            break;
        case 'g': /* gpu mode */
            break;
        default:
            break;
    }
    return 0;
}
static int Preamble(Control *control)
{
    fprintf(stdout, "**********************************************************\n");
    fprintf(stdout, "*                        ArtraCFD                        *\n");
    fprintf(stdout, "*                    <By Huangrui Mo>                    *\n");
    fprintf(stdout, "**********************************************************\n\n");
    fprintf(stdout, "Enter 'help' for more information\n");
    fprintf(stdout, "**********************************************************\n\n");
    String currentLine = {'\0'}; /* store the current read line */
    while (1) {
        fprintf(stdout, "\nArtraCFD << ");
        Fgets(currentLine, sizeof currentLine, stdin); /* read a line */
        CommandLineProcessor(currentLine); /* process current line */
        fprintf(stdout, "\n");
        if (0 == strncmp(currentLine, "help", sizeof currentLine)) {
            fprintf(stdout, "Options under interactive environment:\n\n");
            fprintf(stdout, "[help]    show this information\n");
            fprintf(stdout, "[init]    generate files for a sample case\n");
            fprintf(stdout, "[solve]   solve current case in serial mode\n");
            fprintf(stdout, "[calc]    access expression calculator\n");
            fprintf(stdout, "[manual]  show a brief user manual of ArtraCFD\n");
            fprintf(stdout, "[exit]    exit program\n");
            continue;
        }
        if (0 == strncmp(currentLine, "init", sizeof currentLine)) {
            GenerateCaseFiles();
            fprintf(stdout, "a sample case generated successfully\n");
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
            control->runMode = 's'; /* change mode to serial */
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
    fprintf(stdout, "        -m runmode        run mode: interact, serial, threaded, mpi, gpu\n");
    fprintf(stdout, "        -n nprocessors    number of processors\n");
    fprintf(stdout, "NOTES:\n");
    fprintf(stdout, "        default run mode is 'interact'\n");
    return 0;
}
/* a good practice: end file with a newline */

