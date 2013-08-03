
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = Classic                  # Clawpack package to use
CLAW_EXE = xclaw                    # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS =


# ---------------------------------
# List of sources for this program:
# ---------------------------------


# Clawpack library to be used:
CLAWLIB = $(CLAW)/classic/src/1d
    
SOURCES = \
  qinit.f \
  src1.f \
  setprob.f \
  $(CLAW)/riemann/src/rp1_acoustics.f90 \
  $(CLAWLIB)/setaux.f \
  $(CLAWLIB)/claw1ez.f \
  $(CLAWLIB)/b4step1.f \
  $(CLAWLIB)/bc1.f \
  $(CLAWLIB)/out1.f \
  $(CLAWLIB)/claw1.f \
  $(CLAWLIB)/step1.f \
  $(CLAWLIB)/copyq1.f \
  $(CLAWLIB)/limiter.f \
  $(CLAWLIB)/philim.f\
  $(CLAWLIB)/opendatafile.f

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###
    
