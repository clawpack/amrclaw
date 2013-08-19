
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAWUTIL)/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
EXE = xamr                          # Executable to create
# 3D setrun not set up yet
#SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS ?=  -g -fbounds-check -ffpe-trap=zero,overflow,invalid -fopenmp

# ---------------------------------
# List of sources for this program:
# ---------------------------------

MODULES = \
   $(AMRCLAW)/src/3d/amr_module.f90 \

SOURCES = \
  qinit.f \
  setprob.f \
  ../rp/rpn3eu.f \
  ../rp/rpt3eu.f \
  ../rp/rptt3eu.f \
  $(AMRCLAW)/src/3d/bc3amr.f \
  $(AMRCLAW)/src/3d/setaux.f \
  $(AMRCLAW)/src/3d/amr3.f90 \
  $(AMRCLAW)/src/3d/opendatafile.f \
  $(AMRCLAW)/src/3d/b4step3.f \
  $(AMRCLAW)/src/3d/qad.f \
  $(AMRCLAW)/src/3d/src3.f \
  $(AMRCLAW)/src/3d/src1d.f \
  $(AMRCLAW)/src/3d/advanc.f \
  $(AMRCLAW)/src/3d/bound.f \
  $(AMRCLAW)/src/3d/stepgrid.f \
  $(AMRCLAW)/src/3d/auxcoarsen.f \
  $(AMRCLAW)/src/3d/fixcapaq.f \
  $(AMRCLAW)/src/3d/estdt.f \
  $(AMRCLAW)/src/3d/igetsp.f \
  $(AMRCLAW)/src/3d/reclam.f \
  $(AMRCLAW)/src/3d/birect.f \
  $(AMRCLAW)/src/3d/cleanup.f \
  $(AMRCLAW)/src/3d/colate.f \
  $(AMRCLAW)/src/3d/errest.f \
  $(AMRCLAW)/src/3d/flag2refine.f \
  $(AMRCLAW)/src/3d/allowflag.f \
  $(AMRCLAW)/src/3d/bufnst.f \
  $(AMRCLAW)/src/3d/spest.f \
  $(AMRCLAW)/src/3d/errf1.f \
  $(AMRCLAW)/src/3d/gfixup.f \
  $(AMRCLAW)/src/3d/filval.f \
  $(AMRCLAW)/src/3d/filpatch.f \
  $(AMRCLAW)/src/3d/prefilp.f \
  $(AMRCLAW)/src/3d/flglvl.f \
  $(AMRCLAW)/src/3d/fluxad.f \
  $(AMRCLAW)/src/3d/fluxsv.f \
  $(AMRCLAW)/src/3d/ginit.f \
  $(AMRCLAW)/src/3d/grdfit.f \
  $(AMRCLAW)/src/3d/intfil.f \
  $(AMRCLAW)/src/3d/moment.f \
  $(AMRCLAW)/src/3d/nestck.f \
  $(AMRCLAW)/src/3d/prepc.f \
  $(AMRCLAW)/src/3d/prepf.f \
  $(AMRCLAW)/src/3d/projec.f \
  $(AMRCLAW)/src/3d/signs.f \
  $(AMRCLAW)/src/3d/findcut.f \
  $(AMRCLAW)/src/3d/smartbis.f \
  $(AMRCLAW)/src/3d/putnod.f \
  $(AMRCLAW)/src/3d/putsp.f \
  $(AMRCLAW)/src/3d/regrid.f \
  $(AMRCLAW)/src/3d/setgrd.f \
  $(AMRCLAW)/src/3d/setuse.f \
  $(AMRCLAW)/src/3d/stst1.f \
  $(AMRCLAW)/src/3d/tick.f \
  $(AMRCLAW)/src/3d/trimbd.f \
  $(AMRCLAW)/src/3d/update.f \
  $(AMRCLAW)/src/3d/nodget.f \
  $(AMRCLAW)/src/3d/upbnd.f  \
  $(AMRCLAW)/src/3d/basic.f \
  $(AMRCLAW)/src/3d/outval.f \
  $(AMRCLAW)/src/3d/copysol.f \
  $(AMRCLAW)/src/3d/outvar.f \
  $(AMRCLAW)/src/3d/outmsh.f \
  $(AMRCLAW)/src/3d/outtre.f \
  $(AMRCLAW)/src/3d/domain.f  \
  $(AMRCLAW)/src/3d/setflags.f \
  $(AMRCLAW)/src/3d/shiftset.f \
  $(AMRCLAW)/src/3d/conck.f \
  $(AMRCLAW)/src/3d/domshrink.f \
  $(AMRCLAW)/src/3d/domprep.f \
  $(AMRCLAW)/src/3d/domup.f \
  $(AMRCLAW)/src/3d/domcopy.f \
  $(AMRCLAW)/src/3d/coarsen.f \
  $(AMRCLAW)/src/3d/intcopy.f \
  $(AMRCLAW)/src/3d/preintcopy.f \
  $(AMRCLAW)/src/3d/icall.f \
  $(AMRCLAW)/src/3d/preicall.f \
  $(AMRCLAW)/src/3d/step3.f \
  $(AMRCLAW)/src/3d/flux3.f \
  $(AMRCLAW)/src/3d/limiter.f \
  $(AMRCLAW)/src/3d/philim.f \
  $(AMRCLAW)/src/3d/cstore.f \
  $(AMRCLAW)/src/3d/saveqc.f \
  $(AMRCLAW)/src/3d/valout.f \
  $(AMRCLAW)/src/3d/check.f \
  $(AMRCLAW)/src/3d/restrt.f \
  $(AMRCLAW)/src/3d/init_alloc.f90 \
  $(AMRCLAW)/src/3d/resize_alloc.f90 \
  $(AMRCLAW)/src/3d/restrt_alloc.f90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

