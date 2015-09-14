
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
RESTART ?= False                    # Should = clawdata.restart in setrun

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?=  


# ---------------------------------
# List of sources for this program:
# ---------------------------------

AMRLIB = $(CLAW)/amrclaw/src/2d

MODULES = \
  $(AMRLIB)/amr_module.f90 \
  $(AMRLIB)/regions_module.f90 \
  $(AMRLIB)/gauges_module.f90

SOURCES = \
  qinit.f \
  qtrue.f \
  setprob.f \
  bc2amr.f \
  flag2refine2.f90 \
  $(CLAW)/riemann/src/rpn2_advection.f \
  $(CLAW)/riemann/src/rpt2_advection.f \
  $(AMRLIB)/amr2.f90 \
  $(AMRLIB)/setaux.f90 \
  $(AMRLIB)/domain.f  \
  $(AMRLIB)/b4step2.f90 \
  $(AMRLIB)/qad.f \
  $(AMRLIB)/src2.f90 \
  $(AMRLIB)/src1d.f90 \
  $(AMRLIB)/advanc.f \
  $(AMRLIB)/bound.f90 \
  $(AMRLIB)/stepgrid.f \
  $(AMRLIB)/stepgrid_dimSplit.f \
  $(AMRLIB)/auxcoarsen.f \
  $(AMRLIB)/fixcapaq.f \
  $(AMRLIB)/estdt.f \
  $(AMRLIB)/init_iflags.f \
  $(AMRLIB)/igetsp.f \
  $(AMRLIB)/reclam.f \
  $(AMRLIB)/birect.f \
  $(AMRLIB)/cleanup.f \
  $(AMRLIB)/colate2.f \
  $(AMRLIB)/bufnst2.f \
  $(AMRLIB)/spest2.f \
  ./valout.f \
  $(AMRLIB)/flagregions2.f90 \
  $(AMRLIB)/errest.f \
  $(AMRLIB)/errf1.f \
  $(AMRLIB)/gfixup.f \
  $(AMRLIB)/filval.f90 \
  $(AMRLIB)/filpatch.f90 \
  $(AMRLIB)/prefilp.f90 \
  $(AMRLIB)/flglvl2.f \
  $(AMRLIB)/flagger.f \
  $(AMRLIB)/prepregstep.f \
  $(AMRLIB)/prepbigstep.f \
  $(AMRLIB)/fluxad.f \
  $(AMRLIB)/fluxsv.f \
  $(AMRLIB)/ginit.f \
  $(AMRLIB)/grdfit2.f \
  $(AMRLIB)/intfil.f \
  $(AMRLIB)/moment.f \
  $(AMRLIB)/nestck2.f \
  $(AMRLIB)/prepf.f \
  $(AMRLIB)/prepc.f \
  $(AMRLIB)/projec2.f \
  $(AMRLIB)/signs.f \
  $(AMRLIB)/findcut.f \
  $(AMRLIB)/smartbis.f \
  $(AMRLIB)/putnod.f \
  $(AMRLIB)/putsp.f \
  $(AMRLIB)/regrid.f \
  $(AMRLIB)/setgrd.f \
  $(AMRLIB)/setuse.f \
  $(AMRLIB)/stst1.f \
  $(AMRLIB)/tick.f \
  $(AMRLIB)/trimbd.f90 \
  $(AMRLIB)/update.f \
  $(AMRLIB)/nodget.f \
  $(AMRLIB)/upbnd.f  \
  $(AMRLIB)/basic.f \
  $(AMRLIB)/outval.f \
  $(AMRLIB)/copysol.f \
  $(AMRLIB)/outvar.f \
  $(AMRLIB)/outmsh.f \
  $(AMRLIB)/outtre.f \
  $(AMRLIB)/setflags.f \
  $(AMRLIB)/shiftset2.f \
  $(AMRLIB)/conck.f \
  $(AMRLIB)/domshrink.f \
  $(AMRLIB)/domprep.f \
  $(AMRLIB)/domup.f \
  $(AMRLIB)/domcopy.f \
  $(AMRLIB)/coarsen.f \
  $(AMRLIB)/intcopy.f \
  $(AMRLIB)/preintcopy.f \
  $(AMRLIB)/icall.f \
  $(AMRLIB)/preicall.f \
  $(AMRLIB)/step2.f90 \
  $(AMRLIB)/step2x.f90 \
  $(AMRLIB)/step2y.f90 \
  $(AMRLIB)/flux2.f \
  $(AMRLIB)/flux2_dimSplit.f \
  $(AMRLIB)/inlinelimiter.f \
  $(AMRLIB)/cstore.f \
  $(AMRLIB)/saveqc.f \
  $(AMRLIB)/check.f \
  $(AMRLIB)/restrt.f \
  $(AMRLIB)/quick_sort1.f \
  $(AMRLIB)/opendatafile.f \
  $(AMRLIB)/addflags.f \
  $(AMRLIB)/domgrid.f \
  $(AMRLIB)/drivesort.f \
  $(AMRLIB)/flagcheck.f \
  $(AMRLIB)/setdomflags.f \
  $(AMRLIB)/setIndices.f \
  $(AMRLIB)/setPhysBndry.f \
  $(AMRLIB)/setPhysBndryFlags.f \
  $(AMRLIB)/coarseGridFlagSet.f \
  $(AMRLIB)/griddomcopy.f \
  $(AMRLIB)/griddomshrink.f \
  $(AMRLIB)/griddomup.f \
  $(AMRLIB)/baseCheck.f \
  $(AMRLIB)/init_alloc.f90 \
  $(AMRLIB)/restrt_alloc.f90 \
  $(AMRLIB)/resize_alloc.f90 \
# $(AMRLIB)/resize_alloc_static.f90 \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

