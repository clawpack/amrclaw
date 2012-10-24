#
#               Makefile for the amr code:
#
#       To compile the AMR library files type:    make amrclawlib
#       (from the CLAW sub-directory)
#
#       To make an executable, type:              make xamr
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xamrhdf
#       (from the application sub-directory)
#
#       To compile a single file.f type:          make file.o
#       (from the application sub-directory)
#
#       To combine source + include files type:   make backup
#       (from the application sub-directory
#        date is appended)
#
FFLAGS = -g -fbounds-check -ffpe-trap=invalid,overflow,zero 
LFLAGS = -g
F77     = gfortran -c
LINK    = gfortran

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

LIBOBJECTS = \
  $(CLAW)/amrclaw/src/3d/amr3ez.o \
  $(CLAW)/amrclaw/src/3d/bc3amr.o \
  $(CLAW)/amrclaw/src/3d/setaux.o \
  $(CLAW)/amrclaw/src/3d/setprob.o \
  $(CLAW)/amrclaw/src/3d/b4step3.o \
  $(CLAW)/amrclaw/src/3d/advanc.o \
  $(CLAW)/amrclaw/src/3d/bound.o \
  $(CLAW)/amrclaw/src/3d/stepgrid.o \
  $(CLAW)/amrclaw/src/3d/estdt.o \
  $(CLAW)/amrclaw/src/3d/flag2refine.o  \
  $(CLAW)/amrclaw/src/3d/allowflag.o  \
  $(CLAW)/amrclaw/src/3d/spest.o  \
  $(CLAW)/amrclaw/src/3d/igetsp.o \
  $(CLAW)/amrclaw/src/3d/reclam.o \
  $(CLAW)/amrclaw/src/3d/birect.o \
  $(CLAW)/amrclaw/src/3d/cleanup.o \
  $(CLAW)/amrclaw/src/3d/colate.o \
  $(CLAW)/amrclaw/src/3d/errest.o \
  $(CLAW)/amrclaw/src/3d/bufnst.o \
  $(CLAW)/amrclaw/src/3d/errf1.o \
  $(CLAW)/amrclaw/src/3d/gfixup.o \
  $(CLAW)/amrclaw/src/3d/filval.o \
  $(CLAW)/amrclaw/src/3d/fixcapaq.o \
  $(CLAW)/amrclaw/src/3d/filpatch.o \
  $(CLAW)/amrclaw/src/3d/filpatch2.o \
  $(CLAW)/amrclaw/src/3d/filpatch3.o \
  $(CLAW)/amrclaw/src/3d/prefilp.o \
  $(CLAW)/amrclaw/src/3d/prefil2.o \
  $(CLAW)/amrclaw/src/3d/prefil3.o \
  $(CLAW)/amrclaw/src/3d/flglvl.o \
  $(CLAW)/amrclaw/src/3d/fluxad.o \
  $(CLAW)/amrclaw/src/3d/fluxsv.o \
  $(CLAW)/amrclaw/src/3d/ginit.o \
  $(CLAW)/amrclaw/src/3d/grdfit.o \
  $(CLAW)/amrclaw/src/3d/intfil.o \
  $(CLAW)/amrclaw/src/3d/moment.o \
  $(CLAW)/amrclaw/src/3d/nestck.o \
  $(CLAW)/amrclaw/src/3d/prepc.o \
  $(CLAW)/amrclaw/src/3d/prepf.o \
  $(CLAW)/amrclaw/src/3d/projec.o \
  $(CLAW)/amrclaw/src/3d/signs.o \
  $(CLAW)/amrclaw/src/3d/findcut.o \
  $(CLAW)/amrclaw/src/3d/smartbis.o \
  $(CLAW)/amrclaw/src/3d/putnod.o \
  $(CLAW)/amrclaw/src/3d/putsp.o \
  $(CLAW)/amrclaw/src/3d/regrid.o \
  $(CLAW)/amrclaw/src/3d/setgrd.o \
  $(CLAW)/amrclaw/src/3d/setuse.o \
  $(CLAW)/amrclaw/src/3d/stst1.o \
  $(CLAW)/amrclaw/src/3d/tick.o \
  $(CLAW)/amrclaw/src/3d/trimbd.o \
  $(CLAW)/amrclaw/src/3d/update.o \
  $(CLAW)/amrclaw/src/3d/nodget.o \
  $(CLAW)/amrclaw/src/3d/upbnd.o  \
  $(CLAW)/amrclaw/src/3d/basic.o \
  $(CLAW)/amrclaw/src/3d/outval.o \
  $(CLAW)/amrclaw/src/3d/copysol.o \
  $(CLAW)/amrclaw/src/3d/outvar.o \
  $(CLAW)/amrclaw/src/3d/outmsh.o \
  $(CLAW)/amrclaw/src/3d/outtre.o \
  $(CLAW)/amrclaw/src/3d/domain.o  \
  $(CLAW)/amrclaw/src/3d/setflags.o \
  $(CLAW)/amrclaw/src/3d/shiftset.o \
  $(CLAW)/amrclaw/src/3d/conck.o \
  $(CLAW)/amrclaw/src/3d/domshrink.o \
  $(CLAW)/amrclaw/src/3d/domprep.o \
  $(CLAW)/amrclaw/src/3d/domup.o \
  $(CLAW)/amrclaw/src/3d/domcopy.o \
  $(CLAW)/amrclaw/src/3d/coarsen.o \
  $(CLAW)/amrclaw/src/3d/auxcoarsen.o \
  $(CLAW)/amrclaw/src/3d/intcopy.o \
  $(CLAW)/amrclaw/src/3d/preintcopy.o \
  $(CLAW)/amrclaw/src/3d/icall.o \
  $(CLAW)/amrclaw/src/3d/preicall.o \
  $(CLAW)/amrclaw/src/3d/flux3.o \
  $(CLAW)/amrclaw/src/3d/flux3fw.o \
  $(CLAW)/amrclaw/src/3d/inlinelimiter.o \
  $(CLAW)/amrclaw/src/3d/limiter.o \
  $(CLAW)/amrclaw/src/3d/philim.o \
  $(CLAW)/amrclaw/src/3d/cstore.o \
  $(CLAW)/amrclaw/src/3d/qad.o \
  $(CLAW)/amrclaw/src/3d/src3.o \
  $(CLAW)/amrclaw/src/3d/src1d.o \
  $(CLAW)/amrclaw/src/3d/saveqc.o \
  $(CLAW)/amrclaw/src/3d/step3.o \
  $(CLAW)/amrclaw/src/3d/check.o \
  $(CLAW)/amrclaw/src/3d/restrt.o \
  $(CLAW)/amrclaw/src/3d/valout.o \
  $(CLAW)/amrclaw/src/3d/check_hdf.o \
  $(CLAW)/amrclaw/src/3d/restrt_hdf.o \
  $(CLAW)/amrclaw/src/3d/valout_hdf.o

SOURCES = \
  $(CLAW)/amrclaw/src/3d/amr3ez.f \
  $(CLAW)/amrclaw/src/3d/bc3amr.f \
  $(CLAW)/amrclaw/src/3d/setaux.f \
  $(CLAW)/amrclaw/src/3d/setprob.f \
  $(CLAW)/amrclaw/src/3d/b4step3.f \
  $(CLAW)/amrclaw/src/3d/advanc.f \
  $(CLAW)/amrclaw/src/3d/bound.f \
  $(CLAW)/amrclaw/src/3d/stepgrid.f \
  $(CLAW)/amrclaw/src/3d/estdt.f \
  $(CLAW)/amrclaw/src/3d/flag2refine.f  \
  $(CLAW)/amrclaw/src/3d/allowflag.f  \
  $(CLAW)/amrclaw/src/3d/spest.f  \
  $(CLAW)/amrclaw/src/3d/igetsp.f \
  $(CLAW)/amrclaw/src/3d/reclam.f \
  $(CLAW)/amrclaw/src/3d/birect.f \
  $(CLAW)/amrclaw/src/3d/cleanup.f \
  $(CLAW)/amrclaw/src/3d/colate.f \
  $(CLAW)/amrclaw/src/3d/bufnst.f \
  $(CLAW)/amrclaw/src/3d/errest.f \
  $(CLAW)/amrclaw/src/3d/errf1.f \
  $(CLAW)/amrclaw/src/3d/gfixup.f \
  $(CLAW)/amrclaw/src/3d/filval.f \
  $(CLAW)/amrclaw/src/3d/fixcapaq.f \
  $(CLAW)/amrclaw/src/3d/filpatch.f \
  $(CLAW)/amrclaw/src/3d/filpatch2.f \
  $(CLAW)/amrclaw/src/3d/filpatch3.f \
  $(CLAW)/amrclaw/src/3d/prefilp.f \
  $(CLAW)/amrclaw/src/3d/prefil2.f \
  $(CLAW)/amrclaw/src/3d/prefil3.f \
  $(CLAW)/amrclaw/src/3d/flglvl.f \
  $(CLAW)/amrclaw/src/3d/fluxad.f \
  $(CLAW)/amrclaw/src/3d/fluxsv.f \
  $(CLAW)/amrclaw/src/3d/ginit.f \
  $(CLAW)/amrclaw/src/3d/grdfit.f \
  $(CLAW)/amrclaw/src/3d/intfil.f \
  $(CLAW)/amrclaw/src/3d/moment.f \
  $(CLAW)/amrclaw/src/3d/nestck.f \
  $(CLAW)/amrclaw/src/3d/prepc.f \
  $(CLAW)/amrclaw/src/3d/prepf.f \
  $(CLAW)/amrclaw/src/3d/projec.f \
  $(CLAW)/amrclaw/src/3d/signs.f \
  $(CLAW)/amrclaw/src/3d/findcut.f \
  $(CLAW)/amrclaw/src/3d/smartbis.f \
  $(CLAW)/amrclaw/src/3d/putnod.f \
  $(CLAW)/amrclaw/src/3d/putsp.f \
  $(CLAW)/amrclaw/src/3d/regrid.f \
  $(CLAW)/amrclaw/src/3d/setgrd.f \
  $(CLAW)/amrclaw/src/3d/setuse.f \
  $(CLAW)/amrclaw/src/3d/stst1.f \
  $(CLAW)/amrclaw/src/3d/tick.f \
  $(CLAW)/amrclaw/src/3d/trimbd.f \
  $(CLAW)/amrclaw/src/3d/update.f \
  $(CLAW)/amrclaw/src/3d/nodget.f \
  $(CLAW)/amrclaw/src/3d/upbnd.f  \
  $(CLAW)/amrclaw/src/3d/basic.f \
  $(CLAW)/amrclaw/src/3d/outval.f \
  $(CLAW)/amrclaw/src/3d/copysol.f \
  $(CLAW)/amrclaw/src/3d/outvar.f \
  $(CLAW)/amrclaw/src/3d/outmsh.f \
  $(CLAW)/amrclaw/src/3d/outtre.f \
  $(CLAW)/amrclaw/src/3d/domain.f  \
  $(CLAW)/amrclaw/src/3d/setflags.f \
  $(CLAW)/amrclaw/src/3d/shiftset.f \
  $(CLAW)/amrclaw/src/3d/conck.f \
  $(CLAW)/amrclaw/src/3d/domshrink.f \
  $(CLAW)/amrclaw/src/3d/domprep.f \
  $(CLAW)/amrclaw/src/3d/domup.f \
  $(CLAW)/amrclaw/src/3d/domcopy.f \
  $(CLAW)/amrclaw/src/3d/coarsen.f \
  $(CLAW)/amrclaw/src/3d/auxcoarsen.f \
  $(CLAW)/amrclaw/src/3d/intcopy.f \
  $(CLAW)/amrclaw/src/3d/preintcopy.f \
  $(CLAW)/amrclaw/src/3d/icall.f \
  $(CLAW)/amrclaw/src/3d/preicall.f \
  $(CLAW)/amrclaw/src/3d/step3.f \
  $(CLAW)/amrclaw/src/3d/flux3.f \
  $(CLAW)/amrclaw/src/3d/flux3fw.f \
  $(CLAW)/amrclaw/src/3d/inlinelimiter.f \
  $(CLAW)/amrclaw/src/3d/limiter.f \
  $(CLAW)/amrclaw/src/3d/philim.f \
  $(CLAW)/amrclaw/src/3d/cstore.f \
  $(CLAW)/amrclaw/src/3d/qad.f \
  $(CLAW)/amrclaw/src/3d/src3.f \
  $(CLAW)/amrclaw/src/3d/src1d.f \
  $(CLAW)/amrclaw/src/3d/saveqc.f \
  $(CLAW)/amrclaw/src/3d/check.f \
  $(CLAW)/amrclaw/src/3d/restrt.f \
  $(CLAW)/amrclaw/src/3d/valout.f \
  $(CLAW)/amrclaw/src/3d/check_hdf.f \
  $(CLAW)/amrclaw/src/3d/restrt_hdf.f \
  $(CLAW)/amrclaw/src/3d/valout_hdf.f

IFS = \
#  cuser.i \
  $(CLAW)/amrclaw/src/3d/call.i

amrclawlib: $(LIBOBJECTS)

clean:
	-rm  -f *.o

new:
				make clean
				make amrclawlib

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > amrclaw.`date +%y%m%d`

program : $(SOURCES)
	cat $(SOURCES) > amrclaw.f

### DO NOT remove this line - make depends on it ###
