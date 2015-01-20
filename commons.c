/****************************************************************************
 * Common Functions                                                         *
 * Last-modified: 19 Jan 2015 10:02:01 AM
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
/*
 * Command line processor
 * Get rid of end of line, and information after #.
 * Get rid of before and after tabs, replace between tabs with a space.
 * Get rid of before and after spaces, retain only one space in words.
 * If no other information exists, the lineCommand turns to a NULL string.
 */
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
        ++scanner;
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
            if (*(scanner - 1) != ' ') { /* because tabs are all replaced by space */
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
    /* check whether ended with a space */
    if ((receiver != lineCommand) && *(receiver - 1) == ' ') {
        *(receiver -1) = '\0';
    }
    *receiver = '\0'; /* add string terminator */
    return 0;
}
/*
 * Fatal error control, print information and then exit.
 * Once the process exits, the operating system is able to 
 * free all dynamically allocated memory associated with the process
 */
void FatalError(const char *statement)
{
    printf("Error: %s!\n", statement);
    printf("\n\n**********************************************************\n\n");
    exit(EXIT_FAILURE); /* indicate failure */
}
/*
 * Show information to terminal
 */
int ShowInformation(const char *statement)
{
    if (strcmp(statement, "Session End") == 0) {
        printf("\n**********************************************************\n\n");
        return 0;
    }
    printf("%s\n", statement);
    return 0;
}
/*
 * Assign storage for a 1-order pointer with any datatype.
 * Note:
 *  - in C, don't need to cast the return value of malloc. The pointer to
 *    void returned by malloc is automatically converted to the correct type.
 *    However, if compile with a C++ compiler, a cast is needed.
 *  - malloc does not initialize the storage; this means that the assigned
 *    memory may contain random or unexpected values!
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
void *AssignStorage(const int idxMax, const char *type)
{
    void *pointer = NULL;
    if (strcmp(type, "double") == 0) {
        pointer = malloc(idxMax * sizeof(double));
        if (pointer == NULL) {
            FatalError("Memory allocation failed");
        }
        return pointer;
    }
    if (strcmp(type, "char") == 0) {
        pointer = malloc(idxMax * sizeof(char));
        if (pointer == NULL) {
            FatalError("Memory allocation failed");
        }
        return pointer;
    }
    if (strcmp(type, "int") == 0) {
        pointer = malloc(idxMax * sizeof(int));
        if (pointer == NULL) {
            FatalError("Memory allocation failed");
        }
        return pointer;
    }
    FatalError("Unknown data type");
    return NULL;
}
/*
 * Retrieve storage from a pointer.
 * Note: the original pointer becomes to be a wild pointer after being freed,
 * be aware of this situation! It's a better practice to set pointer back to NULL
 * after calling free, however, this requires passing the address of the pointer.
 */
int RetrieveStorage(void *pointer)
{
    if (pointer != NULL) {
        free(pointer);
    }
    return 0;
}
/* a good practice: end file with a newline */

