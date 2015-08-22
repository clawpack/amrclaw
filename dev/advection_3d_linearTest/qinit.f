c     =====================================================
      subroutine qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower,
     &     dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Sample scalar equation with linear data to test
c     # that second order method stays exact.
c
      implicit double precision (a-h,o-z)
      dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc)
      dimension aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc)

c     # set concentration profile
c     ---------------------------
      do 20 k=1,mz
         zk = zlower + (k-0.5d0)*dz
         do 20 j=1,my
            yj = ylower + (j-0.5d0)*dy
            do 20 i=1,mx
               xi = xlower + (i-0.5d0)*dx
               q(1,i,j,k) = qtrue(xi,yj,zk,0.d0)
   20    continue

      return
      end
