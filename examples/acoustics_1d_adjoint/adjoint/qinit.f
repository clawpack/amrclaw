c     
c     
c=========================================================
      subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c=========================================================
c     
c     # Set initial conditions for adjoint solution
c     # Pulse in pressure, zero velocity, centered around location of interest
c     
c     
      implicit none

      integer, intent(in) :: meqn, mbc, mx, maux
      double precision, intent(in) :: xlower, dx, aux
      double precision, intent(out) :: q
      dimension q(meqn, 1-mbc:mx+mbc)
      dimension aux(maux, 1-mbc:mx+mbc)

      common /cqinit/ beta
      double precision beta

      integer i
      double precision xcell
c     
c     
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx

         q(1,i) = dexp(-beta*(xcell-1.5d0)**2)
         q(2,i) = 0.d0
 150  continue
c
      return
      end
