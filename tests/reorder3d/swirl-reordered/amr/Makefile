#
#               Makefile for the amrclaw code:
#
#       To make an executable, type:            make xamr
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xamrhdf
#       (from the application sub-directory)
#
#       To compile a single file.f type:          make file.o
#       (from the application sub-directory)
#
#       To compile the AMR library files type:    make lib
#       (from the claw/amrclaw/3d/lib sub-directory)
#
#       To combine all source files type:         make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)
#
#       To combine source + include files type:    make backup
#       (from the application sub-directory
#        date is appended)
#
FFLAGS = -O
LFLAGS =
F77    ?= gfortran
LINK   ?= gfortran
HDFLIBS = -L/usr/local/hdf/lib -lmfhdf -ldf -ljpeg -lz

.f.o: ; $(F77) $(FFLAGS) -c $*.f -o $*.o

OBJECTS = \
  ../qinit.o \
  ../../rp-reordered/rpn3adv.o \
  ../../rp-reordered/rpt3adv.o \
  ../../rp-reordered/rptt3adv.o \
  ../setaux.o

LIBOBJECTS = \
  $(CLAW)/amrclaw/3d/lib/amr3ez.o \
  $(CLAW)/amrclaw/3d/lib/bc3amr.o \
  $(CLAW)/amrclaw/3d/lib/b4step3.o \
  $(CLAW)/amrclaw/3d/lib/qad.o \
  $(CLAW)/amrclaw/3d/lib/src3.o \
  $(CLAW)/amrclaw/3d/lib/setprob.o \
  $(CLAW)/amrclaw/3d/lib/src1d.o \
  $(CLAW)/amrclaw/3d/lib/advanc.o \
  $(CLAW)/amrclaw/3d/lib/bound.o \
  $(CLAW)/amrclaw/3d/lib/stepgrid.o \
  $(CLAW)/amrclaw/3d/lib/auxcoarsen.o \
  $(CLAW)/amrclaw/3d/lib/fixcapaq.o \
  $(CLAW)/amrclaw/3d/lib/estdt.o \
  $(CLAW)/amrclaw/3d/lib/igetsp.o \
  $(CLAW)/amrclaw/3d/lib/reclam.o \
  $(CLAW)/amrclaw/3d/lib/birect.o \
  $(CLAW)/amrclaw/3d/lib/cleanup.o \
  $(CLAW)/amrclaw/3d/lib/colate.o \
  $(CLAW)/amrclaw/3d/lib/errest.o \
  $(CLAW)/amrclaw/3d/lib/flag2refine.o \
  $(CLAW)/amrclaw/3d/lib/allowflag.o \
  $(CLAW)/amrclaw/3d/lib/bufnst.o \
  $(CLAW)/amrclaw/3d/lib/spest.o \
  $(CLAW)/amrclaw/3d/lib/errf1.o \
  $(CLAW)/amrclaw/3d/lib/gfixup.o \
  $(CLAW)/amrclaw/3d/lib/filval.o \
  $(CLAW)/amrclaw/3d/lib/filpatch.o \
  $(CLAW)/amrclaw/3d/lib/filpatch2.o \
  $(CLAW)/amrclaw/3d/lib/filpatch3.o \
  $(CLAW)/amrclaw/3d/lib/prefilp.o \
  $(CLAW)/amrclaw/3d/lib/prefil2.o \
  $(CLAW)/amrclaw/3d/lib/prefil3.o \
  $(CLAW)/amrclaw/3d/lib/flglvl.o \
  $(CLAW)/amrclaw/3d/lib/fluxad.o \
  $(CLAW)/amrclaw/3d/lib/fluxsv.o \
  $(CLAW)/amrclaw/3d/lib/ginit.o \
  $(CLAW)/amrclaw/3d/lib/grdfit.o \
  $(CLAW)/amrclaw/3d/lib/intfil.o \
  $(CLAW)/amrclaw/3d/lib/moment.o \
  $(CLAW)/amrclaw/3d/lib/nestck.o \
  $(CLAW)/amrclaw/3d/lib/prepc.o \
  $(CLAW)/amrclaw/3d/lib/prepf.o \
  $(CLAW)/amrclaw/3d/lib/projec.o \
  $(CLAW)/amrclaw/3d/lib/signs.o \
  $(CLAW)/amrclaw/3d/lib/findcut.o \
  $(CLAW)/amrclaw/3d/lib/smartbis.o \
  $(CLAW)/amrclaw/3d/lib/putnod.o \
  $(CLAW)/amrclaw/3d/lib/putsp.o \
  $(CLAW)/amrclaw/3d/lib/regrid.o \
  $(CLAW)/amrclaw/3d/lib/setgrd.o \
  $(CLAW)/amrclaw/3d/lib/setuse.o \
  $(CLAW)/amrclaw/3d/lib/stst1.o \
  $(CLAW)/amrclaw/3d/lib/tick.o \
  $(CLAW)/amrclaw/3d/lib/trimbd.o \
  $(CLAW)/amrclaw/3d/lib/update.o \
  $(CLAW)/amrclaw/3d/lib/nodget.o \
  $(CLAW)/amrclaw/3d/lib/upbnd.o  \
  $(CLAW)/amrclaw/3d/lib/basic.o \
  $(CLAW)/amrclaw/3d/lib/outval.o \
  $(CLAW)/amrclaw/3d/lib/copysol.o \
  $(CLAW)/amrclaw/3d/lib/outvar.o \
  $(CLAW)/amrclaw/3d/lib/outmsh.o \
  $(CLAW)/amrclaw/3d/lib/outtre.o \
  $(CLAW)/amrclaw/3d/lib/domain.o  \
  $(CLAW)/amrclaw/3d/lib/setflags.o \
  $(CLAW)/amrclaw/3d/lib/shiftset.o \
  $(CLAW)/amrclaw/3d/lib/conck.o \
  $(CLAW)/amrclaw/3d/lib/domshrink.o \
  $(CLAW)/amrclaw/3d/lib/domprep.o \
  $(CLAW)/amrclaw/3d/lib/domup.o \
  $(CLAW)/amrclaw/3d/lib/domcopy.o \
  $(CLAW)/amrclaw/3d/lib/coarsen.o \
  $(CLAW)/amrclaw/3d/lib/intcopy.o \
  $(CLAW)/amrclaw/3d/lib/preintcopy.o \
  $(CLAW)/amrclaw/3d/lib/icall.o \
  $(CLAW)/amrclaw/3d/lib/preicall.o \
  $(CLAW)/amrclaw/3d/lib/step3.o \
  $(CLAW)/amrclaw/3d/lib/flux3.o \
  $(CLAW)/amrclaw/3d/lib/limiter.o \
  $(CLAW)/amrclaw/3d/lib/philim.o \
  $(CLAW)/amrclaw/3d/lib/cstore.o \
  $(CLAW)/amrclaw/3d/lib/saveqc.o

SOURCES = \
  ../qinit.f \
  ../../rp-reordered/rpn3adv.f \
  ../../rp-reordered/rpt3adv.f \
  ../../rp-reordered/rptt3adv.f \
  ../setaux.f \
  $(CLAW)/amrclaw/3d/lib/amr3ez.f \
  $(CLAW)/amrclaw/3d/lib/bc3amr.f \
  $(CLAW)/amrclaw/3d/lib/b4step3.f \
  $(CLAW)/amrclaw/3d/lib/qad.f \
  $(CLAW)/amrclaw/3d/lib/src3.f \
  $(CLAW)/amrclaw/3d/lib/setprob.f \
  $(CLAW)/amrclaw/3d/lib/src1d.f \
  $(CLAW)/amrclaw/3d/lib/advanc.f \
  $(CLAW)/amrclaw/3d/lib/bound.f \
  $(CLAW)/amrclaw/3d/lib/stepgrid.f \
  $(CLAW)/amrclaw/3d/lib/auxcoarsen.f \
  $(CLAW)/amrclaw/3d/lib/fixcapaq.f \
  $(CLAW)/amrclaw/3d/lib/estdt.f \
  $(CLAW)/amrclaw/3d/lib/igetsp.f \
  $(CLAW)/amrclaw/3d/lib/reclam.f \
  $(CLAW)/amrclaw/3d/lib/birect.f \
  $(CLAW)/amrclaw/3d/lib/cleanup.f \
  $(CLAW)/amrclaw/3d/lib/colate.f \
  $(CLAW)/amrclaw/3d/lib/errest.f \
  $(CLAW)/amrclaw/3d/lib/flag2refine.f \
  $(CLAW)/amrclaw/3d/lib/allowflag.f \
  $(CLAW)/amrclaw/3d/lib/bufnst.f \
  $(CLAW)/amrclaw/3d/lib/spest.f \
  $(CLAW)/amrclaw/3d/lib/errf1.f \
  $(CLAW)/amrclaw/3d/lib/gfixup.f \
  $(CLAW)/amrclaw/3d/lib/filval.f \
  $(CLAW)/amrclaw/3d/lib/filpatch.f \
  $(CLAW)/amrclaw/3d/lib/filpatch2.f \
  $(CLAW)/amrclaw/3d/lib/filpatch3.f \
  $(CLAW)/amrclaw/3d/lib/prefilp.f \
  $(CLAW)/amrclaw/3d/lib/prefil2.f \
  $(CLAW)/amrclaw/3d/lib/prefil3.f \
  $(CLAW)/amrclaw/3d/lib/flglvl.f \
  $(CLAW)/amrclaw/3d/lib/fluxad.f \
  $(CLAW)/amrclaw/3d/lib/fluxsv.f \
  $(CLAW)/amrclaw/3d/lib/ginit.f \
  $(CLAW)/amrclaw/3d/lib/grdfit.f \
  $(CLAW)/amrclaw/3d/lib/intfil.f \
  $(CLAW)/amrclaw/3d/lib/moment.f \
  $(CLAW)/amrclaw/3d/lib/nestck.f \
  $(CLAW)/amrclaw/3d/lib/prepc.f \
  $(CLAW)/amrclaw/3d/lib/prepf.f \
  $(CLAW)/amrclaw/3d/lib/projec.f \
  $(CLAW)/amrclaw/3d/lib/signs.f \
  $(CLAW)/amrclaw/3d/lib/findcut.f \
  $(CLAW)/amrclaw/3d/lib/smartbis.f \
  $(CLAW)/amrclaw/3d/lib/putnod.f \
  $(CLAW)/amrclaw/3d/lib/putsp.f \
  $(CLAW)/amrclaw/3d/lib/regrid.f \
  $(CLAW)/amrclaw/3d/lib/setgrd.f \
  $(CLAW)/amrclaw/3d/lib/setuse.f \
  $(CLAW)/amrclaw/3d/lib/stst1.f \
  $(CLAW)/amrclaw/3d/lib/tick.f \
  $(CLAW)/amrclaw/3d/lib/trimbd.f \
  $(CLAW)/amrclaw/3d/lib/update.f \
  $(CLAW)/amrclaw/3d/lib/nodget.f \
  $(CLAW)/amrclaw/3d/lib/upbnd.f  \
  $(CLAW)/amrclaw/3d/lib/basic.f \
  $(CLAW)/amrclaw/3d/lib/outval.f \
  $(CLAW)/amrclaw/3d/lib/copysol.f \
  $(CLAW)/amrclaw/3d/lib/outvar.f \
  $(CLAW)/amrclaw/3d/lib/outmsh.f \
  $(CLAW)/amrclaw/3d/lib/outtre.f \
  $(CLAW)/amrclaw/3d/lib/domain.f  \
  $(CLAW)/amrclaw/3d/lib/setflags.f \
  $(CLAW)/amrclaw/3d/lib/shiftset.f \
  $(CLAW)/amrclaw/3d/lib/conck.f \
  $(CLAW)/amrclaw/3d/lib/domshrink.f \
  $(CLAW)/amrclaw/3d/lib/domprep.f \
  $(CLAW)/amrclaw/3d/lib/domup.f \
  $(CLAW)/amrclaw/3d/lib/domcopy.f \
  $(CLAW)/amrclaw/3d/lib/coarsen.f \
  $(CLAW)/amrclaw/3d/lib/intcopy.f \
  $(CLAW)/amrclaw/3d/lib/preintcopy.f \
  $(CLAW)/amrclaw/3d/lib/icall.f \
  $(CLAW)/amrclaw/3d/lib/preicall.f \
  $(CLAW)/amrclaw/3d/lib/step3.f \
  $(CLAW)/amrclaw/3d/lib/flux3.f \
  $(CLAW)/amrclaw/3d/lib/limiter.f \
  $(CLAW)/amrclaw/3d/lib/philim.f \
  $(CLAW)/amrclaw/3d/lib/cstore.f \
  $(CLAW)/amrclaw/3d/lib/saveqc.f

OUTOBJECTS = $(CLAW)/amrclaw/3d/lib/valout.o \
    $(CLAW)/amrclaw/3d/lib/check.o \
    $(CLAW)/amrclaw/3d/lib/restrt.o
OUTSOURCES = $(CLAW)/amrclaw/3d/lib/valout.f \
    $(CLAW)/amrclaw/3d/lib/check.f \
    $(CLAW)/amrclaw/3d/lib/restrt.f

HDFOBJECTS = $(CLAW)/amrclaw/3d/lib/valout_hdf.o \
    $(CLAW)/amrclaw/3d/lib/check_hdf.o \
    $(CLAW)/amrclaw/3d/lib/restrt_hdf.o
HDFSOURCES = $(CLAW)/amrclaw/3d/lib/valout_hdf.f \
    $(CLAW)/amrclaw/3d/lib/check_hdf.f \
    $(CLAW)/amrclaw/3d/lib/restrt_hdf.f

IFS =  $(CLAW)/amrclaw/3d/lib/call.i

xamr: $(OBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xamr

xamrhdf: $(OBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) -o xamr \
	$(HDFLIBS)

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > amrclaw.`date +%y%m%d`

program : $(SOURCES) $(OUTSOURCES)
	cat $(SOURCES) $(OUTSOURCES) > amrclaw.f

programhdf : $(SOURCES) $(HDFSOURCES)
	cat $(SOURCES) $(HDFSOURCES) > amrclaw.f

clean:
	-rm -f $(OBJECTS) xamr xamrhdf

### DO NOT remove this line - make depends on it ###
