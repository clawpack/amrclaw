      subroutine setprob
      implicit double precision (a-h,o-z)
      common /comaux/ z1,c1,z2,c2
c
c     # Set the material parameters for the acoustic equations
c
      open(unit=7,file='setprob.data',status='old',form='formatted')
c
c     # Piecewise constant medium with single interface as specified
c     # in setaux.f

c     # Impedance and sound speed in the two materials:

      read(7,*) z1
      read(7,*) c1
      read(7,*) z2
      read(7,*) c2

      return
      end

