/****************************************************************************
 * Common Functions                                                         *
 * Programmer: Huangrui Mo                                                  *
 * - Follow the Google's C/C++ style Guide.                                 *
 * - This file defines some common operations.                              *
 ****************************************************************************/
/****************************************************************************
 * Required Header Files
 ****************************************************************************/
#include "commons.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <math.h> /* common mathematical functions */
#include <string.h> /* manipulating strings */
/****************************************************************************
 * General functions
 ****************************************************************************/
int CommandLineProcessor(char *lineCommand)
{
    char *scanner = lineCommand; /* copy the address to scanner */
    char *receiver = lineCommand; /* copy the address to receiver */
    /* check whether is a NULL command */
    if (lineCommand[0] == '\0') {
        return 0;
    }
    /* then get rid of before tabs and spaces */
    while ((*scanner == ' ') || (*scanner == '\t')) {
        *scanner = ' '; /* replace tab with space */
        ++scanner; /* update scanner */
    }
    while (*scanner != '\0') { /* until reach the end of original command */
        if ((*scanner == '\r') || (*scanner == '\n')) {
            ++scanner; /* scan the next character */
            continue; /* go to next loop */
        }
        if (*scanner == '#') {
            break; /* no more scan */
        }
        if ((*scanner == ' ') || (*scanner == '\t')) { /* a space or tab */
            if (*(scanner - 1) != ' ') { /* check space, tabs are all replaced */
                /* now its a first space or tab between words */
                *receiver = ' '; /* receive a space */
                ++receiver; /* update the receiver address */
                ++scanner; /* scan the next character */
                continue;
            }
            /* otherwise, do not receive */
            *scanner = ' '; /* replace tab with space */
            ++scanner; /* scan the next character */
            continue;
        }
        /* now its a normal character */
        *receiver = *scanner; /* receive the value of current scanner */
        ++receiver; /* update the receiver address */
        ++scanner; /* scan the next character */
    }
    /* if no useful information in command, receiver didn't receive anything */
    if (receiver == lineCommand) {
        *receiver = '\0';
        return 0;
    }
    /* if received something and also ended with a space */
    if (*(receiver - 1) == ' ') {
        *(receiver -1) = '\0';
        return 0;
    }
    *receiver = '\0'; /* add string terminator */
    return 0;
}
void FatalError(const char *statement)
{
    fprintf(stderr, "error: %s\n", statement);
    exit(EXIT_FAILURE); /* indicate failure */
}
int ShowInformation(const char *statement)
{
    if (strcmp(statement, "Session End") == 0) {
        fprintf(stdout, "\n**********************************************************\n\n");
        return 0;
    }
    fprintf(stdout, "%s\n", statement);
    return 0;
}
/*
 * Assign linear array storage to a pointer with a specific datatype.
 * Note:
 *  - in C, don't need to cast the return value of malloc. The pointer to
 *    void returned by malloc is automatically converted to the correct type.
 *    However, if compile with a C++ compiler, a cast is needed.
 *  - The sizeof operator is used to determine the amount of space a designated
 *    datatype would occupy in memory. To use sizeof, the keyword "sizeof" is
 *    followed by a type name or an expression (which may be merely a variable
 *    name). If a type name is used, it must always be enclosed in parentheses,
 *    whereas expressions can be specified with or without parentheses. 
 *  - The sizeof is a unary operator (not a function!), sizeof gives the size
 *    in units of chars.
 *  - When sizeof is applied to the name of a static array (not allocated
 *    through malloc), the result is the size in bytes (in unit of chars) of the 
 *    whole array, that is, number of elements times the sizeof an array element.
 *    This is one of the few exceptions to the rule that the name of an array is
 *    converted to a pointer to the first element of the array, and is possible
 *    just because the actual array size is fixed and known at compile time, 
 *    when sizeof operator is evaluated.
 *  - When returning a pointer from a function, do not return a pointer that
 *    points to a value that is local to the function or that is a pointer
 *    to a function argument. Pointers to local variables become invalid
 *    when the function exits. In a function, the value returned points to
 *    a static variable or returning a pointer to dynamically allocated 
 *    memory can both be valid.
 */
void *AssignStorage(const int idxMax, const char *dataType)
{
    void *pointer = NULL;
    if (strcmp(dataType, "double") == 0) {
        pointer = malloc(idxMax * sizeof(double));
    }
    if (strcmp(dataType, "Real") == 0) {
        pointer = malloc(idxMax * sizeof(Real));
    }
    if (strcmp(dataType, "float") == 0) {
        pointer = malloc(idxMax * sizeof(float));
    }
    if (strcmp(dataType, "char") == 0) {
        pointer = malloc(idxMax * sizeof(char));
    }
    if (strcmp(dataType, "int") == 0) {
        pointer = malloc(idxMax * sizeof(int));
    }
    if (pointer == NULL) {
        FatalError("memory allocation failed");
    }
    return pointer;
}
int RetrieveStorage(void *pointer)
{
    if (pointer != NULL) {
        free(pointer);
    }
    return 0;
}
/* a good practice: end file with a newline */

