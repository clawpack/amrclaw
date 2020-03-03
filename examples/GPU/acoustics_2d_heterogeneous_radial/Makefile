
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

USE_CUDA = TRUE
FC = pgi
CUDA_PATH ?= /usr/local/cuda-10.0
CUDA_VERSION ?= 10.0
CUDA_CC ?= cc60
 
 
DEBUG = FALSE

CC = pgi
USE_CPP = TRUE

PROFILE = TRUE

USE_OPENMP = TRUE

VERBOSE = FALSE

USE_DOUBLE_PRECISION = TRUE

DIM = 2
NUM_EQUATIONS = 3
NUM_WAVES = 2
NUM_COEFFS = 2
USE_FWAVES = FALSE
USE_CAPA = FALSE

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
  setprob_module.f90

SOURCES = \
  qinit.f \
  setaux.f90 \
  $(CLAW)/riemann/src/rp2_acoustics_variable.f90

ifeq ($(USE_CUDA), TRUE) 
    CPP_SOURCES += stepgrid_cuda.cu
endif
#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

