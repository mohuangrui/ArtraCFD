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
void *AssignStorage(size_t size)
{
    void *pointer = malloc(size);
    if (NULL == pointer) {
        FatalError("memory allocation failed");
    }
    memset(pointer, 0, size); /* initialize to zero */
    return pointer;
}
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

