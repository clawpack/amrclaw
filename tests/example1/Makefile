
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
SETRUN_FILE = setrun.py             # File containing function to make data
OUTDIR = _output                    # Directory for output
SETPLOT_FILE = setplot.py           # File containing function to set plots
PLOTDIR = _plots                    # Directory for plots

OVERWRITE ?= True                   # False ==> make a copy of OUTDIR first
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# intel ifort flags
# FFLAGS =  -g -C -CB -CU -fpe0 -ftrapuv
# Add any desired compiler flags such as -g here:
FFLAGS ?=  -O2
#FFLAGS = -W -Wall -fbounds-check -Wunderflow -O -ffpe-trap=zero,overflow,underflow -g

# ---------------------------------
# List of sources for this program:
# ---------------------------------
MODULES = \
  $(AMRCLAW)/src/2d/amr_module.f90 \

SOURCES = \
  qinit.f \
  setprob.f \
  $(RIEMANN)/src/rpn2_advection.f \
  $(RIEMANN)/src/rpt2_advection.f \
<<<<<<< HEAD
  $(CLAW_LIB)/amr2ez.f \
  $(CLAW_LIB)/bc2amr.f \
  $(CLAW_LIB)/setaux.f \
  $(CLAW_LIB)/b4step2.f \
  $(CLAW_LIB)/qad.f \
  $(CLAW_LIB)/src2.f \
  $(CLAW_LIB)/src1d.f \
  $(CLAW_LIB)/advanc.f \
  $(CLAW_LIB)/bound.f \
  $(CLAW_LIB)/stepgrid.f \
  $(CLAW_LIB)/auxcoarsen.f \
  $(CLAW_LIB)/fixcapaq.f \
  $(CLAW_LIB)/estdt.f \
  $(CLAW_LIB)/igetsp.f \
  $(CLAW_LIB)/reclam.f \
  $(CLAW_LIB)/birect.f \
  $(CLAW_LIB)/cleanup.f \
  $(CLAW_LIB)/colate.f \
  $(CLAW_LIB)/colate2.f \
  $(CLAW_LIB)/bufnst.f \
  $(CLAW_LIB)/bufnst2.f \
  $(CLAW_LIB)/domgrid.f \
  $(CLAW_LIB)/spest2.f \
  $(CLAW_LIB)/flag2refine.f \
  $(CLAW_LIB)/flag2refine2.f \
  $(CLAW_LIB)/allowflag.f \
  $(CLAW_LIB)/errest.f \
  $(CLAW_LIB)/errf1.f \
  $(CLAW_LIB)/gfixup.f \
  $(CLAW_LIB)/filval.f \
  $(CLAW_LIB)/filpatch.f \
  $(CLAW_LIB)/prefilp.f \
  $(CLAW_LIB)/flglvl.f \
  $(CLAW_LIB)/fluxad.f \
  $(CLAW_LIB)/fluxsv.f \
  $(CLAW_LIB)/ginit.f \
  $(CLAW_LIB)/grdfit.f \
  $(CLAW_LIB)/intfil.f \
  $(CLAW_LIB)/moment.f \
  $(CLAW_LIB)/nestck.f \
  $(CLAW_LIB)/nestck2.f \
  $(CLAW_LIB)/flagcheck.f \
  $(CLAW_LIB)/baseCheck.f \
  $(CLAW_LIB)/prepf.f \
  $(CLAW_LIB)/prepc.f \
  $(CLAW_LIB)/projec.f \
  $(CLAW_LIB)/projec2.f \
  $(CLAW_LIB)/signs.f \
  $(CLAW_LIB)/findcut.f \
  $(CLAW_LIB)/smartbis.f \
  $(CLAW_LIB)/putnod.f \
  $(CLAW_LIB)/putsp.f \
  $(CLAW_LIB)/regrid.f \
  $(CLAW_LIB)/setgrd.f \
  $(CLAW_LIB)/setuse.f \
  $(CLAW_LIB)/stst1.f \
  $(CLAW_LIB)/tick.f \
  $(CLAW_LIB)/trimbd.f \
  $(CLAW_LIB)/update.f \
  $(CLAW_LIB)/nodget.f \
  $(CLAW_LIB)/upbnd.f  \
  $(CLAW_LIB)/basic.f \
  $(CLAW_LIB)/outval.f \
  $(CLAW_LIB)/copysol.f \
  $(CLAW_LIB)/outvar.f \
  $(CLAW_LIB)/outmsh.f \
  $(CLAW_LIB)/outtre.f \
  $(CLAW_LIB)/domain.f  \
  $(CLAW_LIB)/setflags.f \
  $(CLAW_LIB)/setdomflags.f \
  $(CLAW_LIB)/setIndices.f \
  $(CLAW_LIB)/coarseGridFlagSet.f \
  $(CLAW_LIB)/addflags.f \
  $(CLAW_LIB)/shiftset.f \
  $(CLAW_LIB)/shiftset2.f \
  $(CLAW_LIB)/conck.f \
  $(CLAW_LIB)/domshrink.f \
  $(CLAW_LIB)/griddomshrink.f \
  $(CLAW_LIB)/domprep.f \
  $(CLAW_LIB)/domup.f \
  $(CLAW_LIB)/griddomup.f \
  $(CLAW_LIB)/domcopy.f \
  $(CLAW_LIB)/griddomcopy.f \
  $(CLAW_LIB)/coarsen.f \
  $(CLAW_LIB)/intcopy.f \
  $(CLAW_LIB)/preintcopy.f \
  $(CLAW_LIB)/icall.f \
  $(CLAW_LIB)/preicall.f \
  $(CLAW_LIB)/step2.f90 \
  $(CLAW_LIB)/flux2.f \
  $(CLAW_LIB)/inlinelimiter.f \
  $(CLAW_LIB)/cstore.f \
  $(CLAW_LIB)/saveqc.f \
  $(CLAW_LIB)/valout.f \
  $(CLAW_LIB)/check.f \
  $(CLAW_LIB)/restrt.f \
  $(CLAW_LIB)/setgauges.f \
  $(CLAW_LIB)/setPhysBndryFlags.f \
  $(CLAW_LIB)/dumpgauge.f \
  $(CLAW_LIB)/quick_sort1.f \
  $(CLAW_LIB)/drivesort.f \
  $(CLAW_LIB)/opendatafile.f \
  $(CLAW_LIB)/init_alloc.f90 \
  $(CLAW_LIB)/restrt_alloc.f90 \
  $(CLAW_LIB)/resize_alloc.f90 \
# $(CLAW_LIB)/resize_alloc_static.f90 \
=======
  $(AMRCLAW)/src/2d/amr2.f \
  $(AMRCLAW)/src/2d/setaux.f \
  $(AMRCLAW)/src/2d/bc2amr.f \
  $(AMRCLAW)/src/2d/b4step2.f \
  $(AMRCLAW)/src/2d/qad.f \
  $(AMRCLAW)/src/2d/src2.f \
  $(AMRCLAW)/src/2d/src1d.f \
  $(AMRCLAW)/src/2d/advanc.f \
  $(AMRCLAW)/src/2d/bound.f \
  $(AMRCLAW)/src/2d/stepgrid.f \
  $(AMRCLAW)/src/2d/auxcoarsen.f \
  $(AMRCLAW)/src/2d/fixcapaq.f \
  $(AMRCLAW)/src/2d/estdt.f \
  $(AMRCLAW)/src/2d/init_iflags.f \
  $(AMRCLAW)/src/2d/igetsp.f \
  $(AMRCLAW)/src/2d/reclam.f \
  $(AMRCLAW)/src/2d/birect.f \
  $(AMRCLAW)/src/2d/cleanup.f \
  $(AMRCLAW)/src/2d/colate.f \
  $(AMRCLAW)/src/2d/bufnst.f \
  $(AMRCLAW)/src/2d/spest.f \
  $(AMRCLAW)/src/2d/flag2refine.f \
  $(AMRCLAW)/src/2d/allowflag.f \
  $(AMRCLAW)/src/2d/errest.f \
  $(AMRCLAW)/src/2d/errf1.f \
  $(AMRCLAW)/src/2d/gfixup.f \
  $(AMRCLAW)/src/2d/filval.f \
  $(AMRCLAW)/src/2d/filpatch.f \
  $(AMRCLAW)/src/2d/prefilp.f \
  $(AMRCLAW)/src/2d/flglvl.f \
  $(AMRCLAW)/src/2d/fluxad.f \
  $(AMRCLAW)/src/2d/fluxsv.f \
  $(AMRCLAW)/src/2d/ginit.f \
  $(AMRCLAW)/src/2d/grdfit.f \
  $(AMRCLAW)/src/2d/intfil.f \
  $(AMRCLAW)/src/2d/moment.f \
  $(AMRCLAW)/src/2d/nestck.f \
  $(AMRCLAW)/src/2d/prepf.f \
  $(AMRCLAW)/src/2d/prepc.f \
  $(AMRCLAW)/src/2d/projec.f \
  $(AMRCLAW)/src/2d/signs.f \
  $(AMRCLAW)/src/2d/findcut.f \
  $(AMRCLAW)/src/2d/smartbis.f \
  $(AMRCLAW)/src/2d/putnod.f \
  $(AMRCLAW)/src/2d/putsp.f \
  $(AMRCLAW)/src/2d/regrid.f \
  $(AMRCLAW)/src/2d/setgrd.f \
  $(AMRCLAW)/src/2d/setuse.f \
  $(AMRCLAW)/src/2d/stst1.f \
  $(AMRCLAW)/src/2d/tick.f \
  $(AMRCLAW)/src/2d/trimbd.f \
  $(AMRCLAW)/src/2d/update.f \
  $(AMRCLAW)/src/2d/nodget.f \
  $(AMRCLAW)/src/2d/upbnd.f  \
  $(AMRCLAW)/src/2d/basic.f \
  $(AMRCLAW)/src/2d/outval.f \
  $(AMRCLAW)/src/2d/copysol.f \
  $(AMRCLAW)/src/2d/outvar.f \
  $(AMRCLAW)/src/2d/outmsh.f \
  $(AMRCLAW)/src/2d/outtre.f \
  $(AMRCLAW)/src/2d/domain.f  \
  $(AMRCLAW)/src/2d/setflags.f \
  $(AMRCLAW)/src/2d/shiftset.f \
  $(AMRCLAW)/src/2d/conck.f \
  $(AMRCLAW)/src/2d/domshrink.f \
  $(AMRCLAW)/src/2d/domprep.f \
  $(AMRCLAW)/src/2d/domup.f \
  $(AMRCLAW)/src/2d/domcopy.f \
  $(AMRCLAW)/src/2d/coarsen.f \
  $(AMRCLAW)/src/2d/intcopy.f \
  $(AMRCLAW)/src/2d/preintcopy.f \
  $(AMRCLAW)/src/2d/icall.f \
  $(AMRCLAW)/src/2d/preicall.f \
  $(AMRCLAW)/src/2d/step2.f90 \
  $(AMRCLAW)/src/2d/flux2.f \
  $(AMRCLAW)/src/2d/inlinelimiter.f \
  $(AMRCLAW)/src/2d/cstore.f \
  $(AMRCLAW)/src/2d/saveqc.f \
  $(AMRCLAW)/src/2d/valout.f \
  $(AMRCLAW)/src/2d/check.f \
  $(AMRCLAW)/src/2d/restrt.f \
  $(AMRCLAW)/src/2d/setgauges.f \
  $(AMRCLAW)/src/2d/dumpgauge.f \
  $(AMRCLAW)/src/2d/quick_sort1.f \
  $(AMRCLAW)/src/2d/opendatafile.f \
  $(AMRCLAW)/src/2d/init_alloc.f90 \
  $(AMRCLAW)/src/2d/restrt_alloc.f90 \
  $(AMRCLAW)/src/2d/resize_alloc.f90 \
# $(AMRCLAW)/src/2d/resize_alloc_static.f90 \
>>>>>>> rjl

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

