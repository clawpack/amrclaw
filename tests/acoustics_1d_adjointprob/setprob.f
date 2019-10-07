      subroutine setprob

      implicit none

      common /cqinit/ beta,ic
      integer ic
      double precision beta

      common /comaux/ Zl, cl, Zr, cr
      double precision rhol, Zl, cl, rhor, Zr, cr
c
c     # Set the material parameters for the acoustic equations
c
      character*25 fname
      character*7 adjointFolder

      fname = 'setprob.data'
      call opendatafile(7, fname)

c     # beta for initial conditions:
      read(7,*) beta
c
c     # Piecewise constant medium with single interface at x=0
c     # Density and sound speed to left and right:
      read(7,*) rhol
      read(7,*) cl
      Zl = rhol*cl

      read(7,*) rhor
      read(7,*) cr
      Zr = rhor*cr

      return
      end
