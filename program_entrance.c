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
static void ConfigureProgram(Control *, Space *);
static void ShowPreamble(Control *);
static void ShowManual(void);
/****************************************************************************
 * Function Definitions
 ****************************************************************************/
int EnterProgram(int argc, char *argv[], Control *control, Space *space)
{
    /*
     * Loop through command line to process options
     * The procedure main takes two arguments: argc and argv. The parameter
     * argc is the number of arguments on the command line (including the
     * program name). The array argv contains the actual arguments.
     * A standard command-line format has the form:
     * command options file1 file2 file3 ...
     * Options are preceded by a dash (-) and are usually a single letter.
     * If the option takes a parameter, it follows the letter with a space.
     * A while loop is used to cycle through the command-line options.
     * One argument always exists: the program name. The expression
     * (argc > 1) checks for additional arguments. The first one is
     * numbered 1. The first character of the first argument is argv[1][0].
     * If this is a dash, there is an option. The switch statement is used
     * to decode the options.
     */
    while ((1 < argc) && ('-' == argv[1][0])) { /* options present */
        if (3 > argc) { /* not enough arguments */
            ShowError("empty entry after: %s\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        switch (argv[1][1]) { /* argv[1][1] is the actual option character */
            /* run mode: -m [gui], [serial], [omp], [mpi], [gpu] */
            case 'm':
                ++argv;
                --argc;
                if (0 == strcmp(argv[1], "gui")) {
                    control->runMode = 'i';
                    break;
                }
                if (0 == strcmp(argv[1], "serial")) {
                    control->runMode = 's';
                    break;
                }
                if (0 == strcmp(argv[1], "omp")) {
                    control->runMode = 'o';
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
                ShowError("bad option: %s\n", argv[1]);
                exit(EXIT_FAILURE);
                /* number of processors: -n nx*ny*nz */
            case 'n':
                ++argv;
                --argc;
                Sscanf(argv[1], 3, "%d*%d*%d", &(control->proc[X]),
                        &(control->proc[Y]), &(control->proc[Z]));
                break;
            default:
                ShowError("bad option: %s\n", argv[1]);
                exit(EXIT_FAILURE);
        }
        /* adjust argument list and count to consume an option */
        ++argv;
        --argc;
    }
    /* check information left */
    if (1 != argc) {
        ShowWarning("unidentified arguments ignored: %s...\n", argv[1]);
    }
    /* configure program according to inputted options */
    ConfigureProgram(control, space);
    return 0;
}
static void ConfigureProgram(Control *control, Space *space)
{
    Partition *const part = &(space->part);
    switch (control->runMode) {
        case 'i': /* gui mode */
            ShowPreamble(control);
            /* fall through */
        case 's': /* serial mode */
            part->proc[X] = 1;
            part->proc[Y] = 1;
            part->proc[Z] = 1;
            part->procN = 1;
            break;
        case 'o': /* omp mode */
            /* fall through */
        case 'm': /* mpi mode */
            part->proc[X] = control->proc[X];
            part->proc[Y] = control->proc[Y];
            part->proc[Z] = control->proc[Z];
            part->procN = control->proc[X] *
                control->proc[Y] * control->proc[Z];
            break;
        case 'g': /* gpu mode */
            break;
        default:
            break;
    }
    return;
}
static void ShowPreamble(Control *control)
{
    ShowInfo("Session");
    ShowInfo("*                         ArtraCFD                         *\n");
    ShowInfo("*                     <By Huangrui Mo>                     *\n");
    ShowInfo("*    Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>     *\n");
    ShowInfo("Session");
    ShowInfo("Enter 'help' for more information\n");
    ShowInfo("Session");
    String str = {'\0'}; /* store the current read line */
    while (1) {
        ShowInfo("\nArtraCFD << ");
        ParseCommand(fgets(str, sizeof str, stdin));
        ShowInfo("\n");
        if (0 == strncmp(str, "help", sizeof str)) {
            ShowInfo("Options in gui environment:\n");
            ShowInfo("[help]    show this information\n");
            ShowInfo("[init]    generate files for a sample case\n");
            ShowInfo("[solve]   solve current case in serial mode\n");
            ShowInfo("[calc]    access expression calculator\n");
            ShowInfo("[manual]  show user manual\n");
            ShowInfo("[exit]    exit program\n");
            continue;
        }
        if (0 == strncmp(str, "init", sizeof str)) {
            GenerateCaseFiles();
            ShowInfo("a sample case generated successfully\n");
            continue;
        }
        if (0 == strncmp(str, "calc", sizeof str)) {
            RunCalculator();
            continue;
        }
        if (0 == strncmp(str, "manual", sizeof str)) {
            ShowManual();
            continue;
        }
        if ('\0' == str[0]) {
            continue;
        }
        if (0 == strncmp(str, "solve", sizeof str)) {
            control->runMode = 's';
            ShowInfo("Session");
            return;
        }
        if (0 == strncmp(str, "exit", sizeof str)) {
            ShowInfo("Session");
            exit(EXIT_SUCCESS);
        }
        /* if non of above is true, then unknow commands */
        ShowWarning("unknown command: %s\n", str);
    }
}
static void ShowManual(void)
{
    ShowInfo("\n            ArtraCFD User Manual\n");
    ShowInfo("SYNOPSIS:\n");
    ShowInfo("        artracfd [-m runmode] [-n nprocessors]\n");
    ShowInfo("OPTIONS:\n");
    ShowInfo("        -m runmode        run mode: gui, serial, omp, mpi, gpu\n");
    ShowInfo("        -n nprocessors    processors per dimension: nx*ny*nz\n");
    ShowInfo("NOTES:\n");
    ShowInfo("        default run mode is gui\n");
    return;
}
/* a good practice: end file with a newline */

