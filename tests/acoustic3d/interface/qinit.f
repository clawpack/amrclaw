
c
c
c
c     =====================================================
      subroutine qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &                 xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
c
      dimension
     & q(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 1-mbc:maxmz+mbc)
      dimension
     & aux(maux, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 1-mbc:maxmz+mbc)
c
c     # set initial data
c     ---------------------------
c
      q = 0.d0

      return
      end
