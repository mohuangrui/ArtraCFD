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
#include <stdarg.h> /* variable-length argument lists */
/****************************************************************************
 * Global Real Constants Definition
 ****************************************************************************/
const Real PI = 3.14159265358979323846;
/****************************************************************************
 * General functions
 ****************************************************************************/
int ParseCommand(char *cmdstr)
{
    char *scanner = cmdstr;
    char *receiver = cmdstr;
    if (NULL == cmdstr) {
        ShowWarning("parse a NULL pointer");
        return 1;
    }
    if ('\0' == *cmdstr) {
        return 1;
    }
    /* skip tabs and spaces at the front */
    while ((' ' == *scanner) || ('\t' == *scanner)) {
        *scanner = ' ';
        ++scanner;
    }
    /* process valid text */
    while ('\0' != *scanner) {
        /* skip special characters */
        if (('\r' == *scanner) || ('\n' == *scanner)) {
            ++scanner;
            continue;
        }
        /* stop at commenting text */
        if ('#' == *scanner) {
            break;
        }
        /* replace tabs with spaces and keep only one space between two words */
        if ((' ' == *scanner) || ('\t' == *scanner)) {
            *scanner = ' ';
            if (' ' != *(scanner - 1)) { /* only check space; tabs are all replaced */
                /* receive the first space between words */
                *receiver = ' ';
                ++receiver;
                ++scanner;
                continue;
            }
            /* otherwise, skip the redundant spaces */
            ++scanner;
            continue;
        }
        /* receive the normal character */
        *receiver = *scanner;
        ++receiver;
        ++scanner;
    }
    /* if received nothing, produce an empty string */
    if (receiver == cmdstr) {
        *receiver = '\0';
        return 1;
    }
    /* if ended with a space, skip it */
    if (' ' == *(receiver - 1)) {
        *(receiver - 1) = '\0';
        return 0;
    }
    *receiver = '\0';
    return 0;
}
char *ParseFormat(char *fmt)
{
    if (sizeof(Real) == sizeof(double)) {
        return fmt;
    }
    char *scanner = fmt;
    char *receiver = fmt;
    while ('\0' != *scanner) {
        if ('l' == *scanner) {
            ++scanner;
            continue;
        }
        *receiver = *scanner;
        ++receiver;
        ++scanner;
    }
    *receiver = '\0';
    return fmt;
}
void ShowError(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    verror("error", fmt, args);
    va_end(args);
    exit(EXIT_FAILURE);
}
void ShowWarning(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    verror("warning", fmt, args);
    va_end(args);
    return;
}
void verror(const char *prefix, const char *fmt, va_list args)
{
    fprintf(stderr, "%s: ", prefix);
    vfprintf(stderr, fmt, args);
    fprintf(stderr, "\n");
    return;
}
void ShowInfo(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    if (0 == strcmp(fmt, "Session")) {
        fprintf(stdout, "************************************************************\n");
    } else {
        vfprintf(stdout, fmt, args);
    }
    va_end(args);
    return;
}
void *AssignStorage(size_t size)
{
    void *pointer = malloc(size);
    if (NULL == pointer) {
        ShowError("memory allocation failed");
    }
    memset(pointer, 0, size);
    return pointer;
}
void RetrieveStorage(void *pointer)
{
    if (NULL != pointer) {
        free(pointer);
    }
    return;
}
void ReadInLine(FILE *fp, const char *line)
{
    String str = {'\0'}; /* store the current read line */
    while (NULL != fgets(str, sizeof str, fp)) {
        ParseCommand(str);
        if (0 == strncmp(str, line, sizeof str)) {
            fseek(fp, 0, SEEK_CUR); /* update offset */
            break;
        }
    }
    return;
}
void WriteToLine(FILE *fp, const char *line)
{
    String str = {'\0'}; /* store the current read line */
    fpos_t pos; /* store position indicator */
    /* find the offset */
    fgetpos(fp, &pos);
    while (NULL != fgets(str, sizeof str, fp)) {
        ParseCommand(str);
        if (0 == strncmp(str, line, sizeof str)) {
            fsetpos(fp, &pos); /* redirect to the target line */
            break;
        }
        fgetpos(fp, &pos);
    }
    return;
}
FILE *Fopen(const char *fname, const char *mode)
{
    FILE *fp = fopen(fname, mode);
    if (NULL == fp) {
        ShowError("failed to open file: %s, mode: %s", fname, mode);
    }
    return fp;
}
void Fread(void *ptr, size_t size, size_t n, FILE *stream)
{
    if (n != fread(ptr, size, n, stream))
    {
        ShowError("fread failed: size: %d, n: %d", size, n);
    }
    return;
}
void Fscanf(FILE *stream, const int n, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    if ((n != vfscanf(stream, fmt, args)) && (0 < n)) {
        ShowError("vfscanf failed: n: %d, fmt: %s", n, fmt);
    }
    va_end(args);
    return;
}
void Sscanf(const char *str, const int n, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    if ((n != vsscanf(str, fmt, args)) && (0 < n)) {
        ShowError("vsscanf failed: str: %s, fmt: %s", str, fmt);
    }
    va_end(args);
    return;
}
void Sread(FILE *stream, const int n, const char *fmt, ...)
{
    String str = {'\0'};
    if (NULL == fgets(str, sizeof str, stream))
    {
        ShowWarning("fgets return a NULL");
    }
    if (0 == n) {
        return;
    }
    va_list args;
    va_start(args, fmt);
    if ((n != vsscanf(str, fmt, args)) && (0 < n)) {
        ShowError("vsscanf failed: str: %s, fmt: %s", str, fmt);
    }
    va_end(args);
    return;
}
/* a good practice: end file with a newline */

