c     
c     
c=========================================================
      subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c=========================================================
c     
c     # Set initial conditions for q.
c     # Pulse in pressure, zero velocity
c     
c     
      implicit none

      integer, intent(in) :: meqn, mbc, mx, maux
      double precision, intent(in) :: xlower, dx, aux
      double precision, intent(out) :: q
      dimension q(meqn, 1-mbc:mx+mbc)
      dimension aux(maux, 1-mbc:mx+mbc)

      common /cqinit/ beta,ic
      double precision beta
      integer ic

      integer i
      double precision xcell
c     
c     
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx

         go to (10,20,30) ic

 10      continue
c     # half ellipse:
         if (xcell.gt.-4d0 .and. xcell.lt.-2d0) then
            q(1,i) = dsqrt(1.d0 - (xcell+3.d0)**2)
         else
            q(1,i) = 0.d0
         endif
         q(2,i) = 0.d0
         go to 150

 20      continue
c     # single discontinuity:
         if (xcell .lt. -2.d0) then
            q(1,i) = 1.d0
         else
            q(1,i) = 0.d0
         endif
         q(2,i) = q(1,i)
         go to 150

 30      continue
c     # Gaussian
         q(1,i) = dexp(-beta*(xcell+2.0d0)**2)  
         q(2,i) = 0.d0
         go to 150

 150  continue
c     
      return
      end
