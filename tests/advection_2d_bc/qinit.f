c     =====================================================
      subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &     dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Sample scalar equation with data that is piecewise constant with
c     # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
c     #     0.1  otherwise
c
      implicit double precision (a-h,o-z)
      dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)

c     # set concentration profile
c     ---------------------------
      do 20 j=1,my
C          yj = ylower + (j-0.5d0)*dy
         do 20 i=1,mx
C             xi = xlower + (i-0.5d0)*dx
C             if (xi.gt.0.1d0 .and. xi.lt.0.6d0 .and.
C      &            yj.gt.0.1d0 .and. yj.lt.0.6d0) then
C                q(1,i,j) = 1.d0
C             else
               q(1,i,j) = 0.1d0
C             endif
   20    continue


      return
      end
