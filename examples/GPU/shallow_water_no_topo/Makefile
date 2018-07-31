
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first

FC = pgi
USE_CUDA = TRUE
# FC = gfortran
# USE_CUDA = FALSE
 
 
# DEBUG = TRUE
DEBUG = FALSE

CC = pgi
USE_CPP = TRUE

PROFILE = TRUE
# PROFILE = FALSE 

USE_OPENMP = TRUE

VERBOSE = FALSE

# USE_DOUBLE_PRECISION = FALSE
USE_DOUBLE_PRECISION = TRUE

DIM = 2
NUM_EQUATIONS = 3
NUM_WAVES = 3
NUM_COEFFS = 0
USE_FWAVES = FALSE

# Environment variable FC should be set to fortran compiler, e.g. gfortran


# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= 

# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d
include $(AMRLIB)/Makefile.amr_2d

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \
  problem_para_module.f90 \
  $(CLAW)/riemann/src/shallow_device_module.f90

SOURCES = \
  qinit.f

ifeq ($(USE_CUDA), TRUE) 
    CPP_SOURCES += stepgrid_cuda.cu
else
endif
#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

