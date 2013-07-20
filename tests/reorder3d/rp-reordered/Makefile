#
#  Makefile for Riemann solvers:
#
#  Type "make" or "make rp" to create .o files needed in
#  applications directories
#
#
FFLAGS =  -c
LFLAGS =
LINK    ?= gfortran
.f.o: ; $(LINK) $(FFLAGS) $*.f

RPOBJECTS = \
  rpn3adv.o \
  rpt3adv.o \
  rptt3adv.o

RPSOURCES = \
  rpn3adv.f \
  rpt3adv.f \
  rptt3adv.f

rp: $(RPOBJECTS)


### DO NOT remove this line - make depends on it ###
