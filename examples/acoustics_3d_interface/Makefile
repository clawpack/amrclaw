
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

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= 


# ---------------------------------
# package sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/3d
include $(AMRLIB)/Makefile.amr_3d

# ---------------------------------------
# package sources specifically to exclude
# (i.e. if a custom replacement source 
#  under a different name is provided)
# Note: There is no need to provide the path to the file; the name alone is sufficient
# ---------------------------------------

EXCLUDE_MODULES = \

EXCLUDE_SOURCES = \
  inlinelimiter.f \

# ----------------------------------------
# List of custom sources for this program:
# ----------------------------------------

MODULES = \

SOURCES = \
  qinit.f \
  setaux.f \
  setprob.f \
  $(CLAW)/riemann/src/rpn3_vc_acoustics.f90 \
  $(CLAW)/riemann/src/rpt3_vc_acoustics.f90 \
  $(CLAW)/riemann/src/rptt3_vc_acoustics.f90 \
  $(AMRLIB)/limiter.f \
  $(AMRLIB)/philim.f \


#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

