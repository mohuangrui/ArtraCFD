#---------------------------------------------------------------------------#
#
#                            ArtraCFD Makefile
#
# Options:
# 'make'             build executable file
# 'make all'         build executable file
# 'make install'     install
# 'make uninstall'   uninstall
# 'make clean'       removes all objects, dependency and executable files
#
# Written by Huangrui Mo, GNU General Public License
#
#---------------------------------------------------------------------------#

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
prefix = ~/MyBin

#
# The directory to install binary executable program
#
bindir = $(prefix)/bin

#
# The directory to install the info files in.
#
infodir = $(prefix)/info

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
#
CC := gcc

#
# Define compiler flags
#   This flag will affect all C compilations uniformly,
#   include implicit rules.
#
#	Shared compiler flags
#    -g        Enable debugging
#    -Wall     Turn on all warnings
#    -Wextra   More restricted warnings
#	GCC compiler flags
#    -O0       No optimization (the default); generates unoptimized code
#              but has the fastest compilation time.
#    -O2       Full optimization; generates highly optimized code and has
#              the slowest compilation time. 
#   ICC compiler flags
#    -O2       Optimize for speed and enable some optimization (default)
#    -O3       Enable all optimizations as O2, and intensive loop optimizations
#    -xP       Enables SSE3, SSE2 and SSE instruction sets optimizations
#
ifeq ($(CC),gcc)
CFLAGS += -Wall -Wextra -O2
else
CFLAGS += -Wall -Wextra -O3 -xP
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
#    Unload intelcc and load gnu module
#
ifeq ($(CC),gcc)
@module unload intel/12.1.3
@module load gcc
else
@module unload gcc
@module load intel/12.1.3
endif

#
# Define any directories containing header files other than /usr/include
#
#    e.g., INCLUDES = -I/home/auxiliary/include  -I./include
#
INCLUDES :=

#
# Define library paths in addition to /usr/lib
#    if wanted libraries not in /usr/lib, specify their path using -Lpath
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
ifeq ($(CC),gcc)
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
