
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

CLAW_PKG = classic
EXE = xclaw
FFLAGS = -g -fbounds-check -ffpe-trap=zero,overflow,invalid
LFLAGS = -g

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


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
  $(CLAW)/classic/src/1d/b4step1.f \
  $(CLAW)/classic/src/1d/setaux.f \
  $(CLAW)/classic/src/1d/bc1.f \
  $(CLAW)/classic/src/1d/claw1ez.f \
  $(CLAW)/classic/src/1d/claw1.f \
  $(CLAW)/classic/src/1d/copyq1.f \
  $(CLAW)/classic/src/1d/driver.f90 \
  $(CLAW)/classic/src/1d/inlinelimiter.f90 \
  $(CLAW)/classic/src/1d/opendatafile.f \
  $(CLAW)/classic/src/1d/out1.f \
  $(CLAW)/classic/src/1d/step1.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)
