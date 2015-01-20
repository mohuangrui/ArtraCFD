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
# Installer
#
INSTALL := install -c
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
#
CC := gcc

#
# Define compiler flags
#   This flag will affect all C compilations uniformly,
#   include implicit rules.
#
#    -g        Enable debugging
#    -Wall     Turn on all warnings
#    -Wextra   More restricted warnings
#    -O        Enable optimization
#
CFLAGS += -Wall -Wextra -O

#
# Preprocessor options
#
CPPFLAGS +=

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
install: all
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

$(DPND): %.d: %.c
	@set -e; rm -f $@; \
		$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
		sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
		rm -f $@.$$$$

# include generated dependencies
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
