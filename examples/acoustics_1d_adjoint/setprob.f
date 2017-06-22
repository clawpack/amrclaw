      subroutine setprob

      use adjoint_module, only: read_adjoint_data, set_time_window
      implicit double precision (a-h,o-z)

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

c     # choice of initial data:
      read(7,*) ic
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

      ! Read data for adjoint problem
      read(7,*) adjointFolder

      ! time period of interest:
      read(7,*) t1
      read(7,*) t2

      call set_time_window(t1, t2)                   !# Set time window
      call read_adjoint_data(trim(adjointFolder))    !# Read adjoint solution

      return
      end
