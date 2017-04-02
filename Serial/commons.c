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
#include "commons.h"
#include <stdio.h> /* standard library for input and output */
#include <stdlib.h> /* dynamic memory allocation and exit */
#include <string.h> /* manipulating strings */
/****************************************************************************
 * General functions
 ****************************************************************************/
int CommandLineProcessor(char *lineCommand)
{
    char *scanner = lineCommand;
    char *receiver = lineCommand;
    /* check whether is a NULL command pointer */
    if (NULL == lineCommand) {
        fprintf(stderr, "warning: process a NULL command line pointer\n");
        return 0;
    }
    /* check whether is a empty command line */
    if ('\0' == *lineCommand) {
        return 0;
    }
    /* then skip tabs and spaces at the front */
    while ((' ' == *scanner) || ('\t' == *scanner)) {
        *scanner = ' '; /* replace tab with space */
        ++scanner; /* update scanner */
    }
    while ('\0' != *scanner) { /* until reach the end of original command */
        if (('\r' == *scanner) || ('\n' == *scanner)) {
            ++scanner; /* scan the next character */
            continue; /* go to next loop */
        }
        if ('#' == *scanner) {
            break; /* no more scan */
        }
        if ((' ' == *scanner) || ('\t' == *scanner)) { /* a space or tab */
            *scanner = ' '; /* replace tab with space */
            if (' ' != *(scanner - 1)) { /* only check space; tabs are all replaced */
                /* now its a first space or tab between words */
                *receiver = ' '; /* receive a space */
                ++receiver; /* update the receiver address */
                ++scanner; /* scan the next character */
                continue;
            }
            /* otherwise, do not receive */
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
    if (' ' == *(receiver - 1)) {
        *(receiver - 1) = '\0';
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
    if (0 == strcmp(statement, "Session End")) {
        fprintf(stdout, "\n**********************************************************\n\n");
        return 0;
    }
    fprintf(stdout, "%s\n", statement);
    return 0;
}
/*
 * Assign linear array storage to a pointer.
 * Note:
 *  - in C, don't need to cast the return value of malloc. The pointer to
 *    void returned by malloc is automatically converted to the correct type.
 *    However, if compile with a C++ compiler, a cast is needed.
 *  - The sizeof operator is used to determine the amount of space a designated
 *    datatype would occupy in memory. To use sizeof, the keyword "sizeof" is
 *    followed by a type name or an expression (which may be merely a variable
 *    name). If a type name is used, it must always be enclosed in parentheses,
 *    whereas expressions can be specified with or without parentheses. 
 *  - The sizeof is a unary operator (not a function), sizeof gives the size
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
 *  - memset treats any memory buffer as a series of bytes, disregarding with
 *    the specific data type. Therefore, it will not set multi-byte types to a
 *    specific non-zero value. For example, int a[100]; memset(a, 1, sizeof(a));
 *    will not set each member of a to the value 1 rather it will set every
 *    byte in the memory buffer taken up by a to 1, which means every four-byte
 *    int is set to the value 0x01010101, which is not the same as 0x00000001.
 */
void *AssignStorage(size_t size)
{
    void *pointer = malloc(size);
    if (NULL == pointer) {
        FatalError("memory allocation failed");
    }
    memset(pointer, 0, size); /* initialize to zero */
    return pointer;
}
/*
 * There are a number of rules you should follow in de-allocating memory.
 * - Prevent access to de-allocated memory. This can be done by setting the
 *   pointer to null after de-allocating, which requires pass the reference
 *   of the target pointer.
 * - De-allocate memory in the reverse order it was allocated. This makes sure
 *   that any dependencies between the allocated memory will not result in 
 *   "dangling pointers". So if one allocated data structure has a pointer to 
 *   another allocated data structure, the second should be de-allocated first.
 * - For a temporary memory block, de-allocate the block before leaving the
 *   routine. If the de-allocation is not done before the routine ends, access
 *   to the memory is lost.
 */
int RetrieveStorage(void *pointer)
{
    if (NULL != pointer) {
        free(pointer);
    }
    return 0;
}
int ReadInLine(FILE *filePointer, const char *lineString)
{
    String currentLine = {'\0'}; /* store the current read line */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, lineString, sizeof currentLine)) {
            break;
        }
    }
    return 0;
}
int WriteToLine(FILE *filePointer, const char *lineString)
{
    String currentLine = {'\0'}; /* store the current read line */
    int offset = 0; /* offset to target line */
    rewind(filePointer); /* seek to the beginning of the file stream */
    while (NULL != fgets(currentLine, sizeof currentLine, filePointer)) {
        CommandLineProcessor(currentLine); /* process current line */
        if (0 == strncmp(currentLine, lineString, sizeof currentLine)) {
            break;
        }
        ++offset;
    }
    /* redirect to the target line */
    rewind(filePointer); /* seek to the beginning of the file stream */
    for (int count = 0; count < offset; ++count) {
        Fgets(currentLine, sizeof currentLine, filePointer);
    }
    return 0;
}
void Fgets(char *str, int num, FILE *stream)
{
    if (NULL == fgets(str, num, stream))
    {
        FatalError("reading information failed...");
    }
    return;
}
void Fread(void *ptr, size_t size, size_t count, FILE *stream)
{
    if (count != fread(ptr, size, count, stream))
    {
        FatalError("reading information failed...");
    }
    return;
}
void VerifyReadConversion(const int num, const int expect)
{
    if (num != expect)
    {
        FatalError("reading information failed...");
    }
    return;
}
/* a good practice: end file with a newline */

