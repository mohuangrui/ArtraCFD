#***************************************************************************#
#                          ArtraCFD Makefile                                #
#                          <By Huangrui Mo>                                 #
# Copyright (C) Huangrui Mo <huangrui.mo@gmail.com>                         #
# This file is part of ArtraCFD.                                            #
# ArtraCFD is free software: you can redistribute it and/or modify it       #
# under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 3 of the License, or         #
# (at your option) any later version.                                       #
#***************************************************************************#

#***************************************************************************#
# Options:
# 'make'             build executable file
# 'make all'         build executable file
# 'make install'     install
# 'make uninstall'   uninstall
# 'make clean'       removes all objects, dependency and executable files
#
# Use 'cat -e -t -v Makefile' to show the presence of tabs with ^I and
# line endings with $, which are vital to ensure that dependencies end
# properly and tabs mark the action for the rules. 
#
#***************************************************************************#

#***************************************************************************#
#
#                           System Configuration
#
#***************************************************************************#

#
# Shell
#
SHELL := /bin/bash

#
# Installer
#
INSTALL := install
INSTALLDATA := $(INSTALL) -m 644

#
# Prefix for each installed program
#
#    e.g., prefix = /usr/local
#
prefix = ~

#
# The directory to install binary executable program
#
bindir = $(prefix)/Bin

#
# The directory to install the info files in.
#
infodir = $(bindir)/info

#
# Define the executable program name
#
BINNAME := artracfd

#
# Path to the source directory, relative to the makefile
#
srcdir = .

#***************************************************************************#
#
#                          Compiler Configuration
#
#***************************************************************************#

#
# Define the compiler
#
#    gcc        GNU C compiler
#    icc        Intel C compiler
#    mpicc      MPI compiler
#
CC := gcc

#
# Define compiler flags
#   This flag affects all C compilations uniformly, include implicit rules.
#
#  Shared compiler flags
#    -Wall     Turn on all warnings
#    -Wextra   More restricted warnings
#    -std=c99 -pedantic  Use ANSI C standard
#    -g        Enable debugging
#    -O0       No optimization; generates unoptimized code for debugging purposes.
#    -O2       Recommended optimization; generates well optimized code. 
#    -O3       Aggressive optimization; recommended for codes that have loops that 
#              heavily use floating-point calculations and process large data sets. 
#  GCC compiler flags
#    -fstrict-aliasing  Allow the compiler to assume the strictest aliasing rules
#              and activates optimizations based on the type of expressions.
#    -Og       Enables optimizations that do not interfere with debugging.
#  ICC compiler flags
#    -ansi-alias  Allow the compiler to assume the strictest aliasing rules 
#              and activates optimizations based on the type of expressions.
#    -fast     Enable processor specific optimization at -O3 level.
#  Use Valgrind for memory access check (http://valgrind.org/)
#    -g -O0    Use this flag to compile the program, then run command line
#    valgrind --leak-check=full --track-origins=yes artracfd arg1 arg2
#  Use Valgrind for cache missing check (http://valgrind.org/)
#    -g -O2    Use this flag to compile the program, then run command line
#    valgrind --tool=cachegrind artracfd arg1 arg2
#  Use google-perftools for performance check (https://github.com/gperftools/gperftools)
#    sudo apt-get update
#    sudo apt-get install google-perftools
#    -g  -O0    Use this flag to compile the program
#    -lprofiler add this flag to the value of LFLAGS
#    make       to generate binary executable, then run command line
#    CPUPROFILE=./cpu.prof ./artracfd arg1 arg2
#    pprof ./artracfd ./cpu.prof       (-pg-like text output)
#    pprof --gv ./artracfd ./cpu.prof  (annotated call-graph via 'gv')
#  Enable floating-point exception handling to debug algorithm
#    trapfpe.c  Add into source code
#    -g  -O0    Use this flag to compile the program
#    gdb ./artracfd
#    (gdb) run -m serial
#    where      Show trace information
#
ifeq ($(CC),icc)
    CFLAGS += -Wall -Wextra -fast -ansi-alias -std=c99 -pedantic
else
    CFLAGS += -Wall -Wextra -O2 -fstrict-aliasing -std=c99 -pedantic
endif

#
# Preprocessor options
#
CPPFLAGS +=

#
# Switch intelcc and gnu module
#
#  When using gcc to compile, it is common to see an error related to
#  <math.h>. This error occurs because the intelcc module is loaded 
#  and is pointing to the intel version of math.h. The Intel version
#  of math.h does not work with the gcc compiler. There are two simple
#  workarounds to fix this problem:
#    Exclusively use icc to compile your jobs.
#    Unload intelcc and load gcc module
#      module unload intelcc
#      module load gcc
#

#
# Define any directories containing header files other than /usr/include
#
#    e.g., INCLUDES = -I/home/auxiliary/include  -I./include
#
INCLUDES :=

#
# Define library paths in addition to /usr/lib
#    If wanted libraries not in /usr/lib, specify their path using -Lpath
#
#    e.g., LFLAGS = -L/home/auxiliary/lib  -L./lib
#
LFLAGS :=

#
# Define any libraries to link into executable, use the -llibname option
#
LIBS := -lm

#***************************************************************************#
#
#                             Make Configuration
#
#***************************************************************************#

#
# Define the C source files
#
SRCS := $(wildcard *.c)

#
# Define the C object files 
#    This uses Suffix Replacement within a macro:
#    $(name:string1=string2)
#    For each word in 'name' replace 'string1' with 'string2'
#    Below we are replacing the suffix .c of all words in the macro SRCS
#    with the .o suffix
#
OBJS := $(SRCS:.c=.o)

#
# Search path for make program
#   make uses VPATH as a search list for both 
#   prerequisites and targets of rules.
#
VPATH :=

#
# Clean list
#
CLEANLIST += $(OBJS) $(BINNAME)

#***************************************************************************#
#
#                                 Build
#
#***************************************************************************#

#
# all
#
.PHONY: all
all: $(BINNAME)
	@echo  $(BINNAME) has been compiled

#
# install
#
.PHONY: install
install:
	@echo "Creating directories"
	@mkdir -p $(bindir)
	@mkdir -p $(infodir)
	@echo "Installing to $(bindir)/$(BINNAME)"
	@$(INSTALL) $(BINNAME) $(bindir)/$(BINNAME)
	@$(INSTALLDATA) $(srcdir)/Makefile $(infodir)

#
# uninstall
#
.PHONY: uninstall
uninstall:
	@echo "Removing  $(bindir)/$(BINNAME)"
	@$(RM)  $(bindir)/$(BINNAME)

#
# Invoke object files
#
$(BINNAME): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) $(CPPFLAGS) -o $@ $(OBJS) $(LFLAGS) $(LIBS)

#
# Static pattern rule for automatic prerequisite generation
#
DPND := $(SRCS:.c=.d)

# Automatic prerequisites flag: -M for any compiler, -MM for GNU to 
# omit system headers. But -MM usually work with ICC without problem.
ifeq ($(CC),icc)
    AUTOPRE := -MM
else
    AUTOPRE := -MM
endif

$(DPND): %.d: %.c
	@set -e; rm -f $@; \
		$(CC) $(AUTOPRE) $(CPPFLAGS) $< > $@.$$$$; \
		sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
		rm -f $@.$$$$

# Include generated dependencies. Extra spaces are allowed and ignored
# at the beginning of the include line, but the first character must 
# NOT be a tab 
ifneq ($(MAKECMDGOALS),clean)
    -include ${DPND}
endif

# add dependency files to clean list
CLEANLIST += $(DPND)

#
# Build object files
#   object files can be built by the automatically generated
#   prerequisites through implicit rules. Generally,
#   to specify additional prerequisites, such as header files,
#   implicit rules are more desirable. 
#

#
# clean
#   When a line starts with ‘@’, the echoing of that line
#   itself is suppressed.
#   a '-' flag makes errors to be ignored
#
.PHONY: clean
clean:
	@echo  cleaning...
	@- $(RM) $(CLEANLIST)

#***************************************************************************#
