c     ============================================
      subroutine b4step3(mbc,mx,my,mz,meqn,q,xlower,ylower,
     &                   zlower,dx,dy,dz,t,dt,maux,aux)
c     ============================================
c
c     # called before each call to step3.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.

c
c     # dummy routine 
c
c     
      implicit double precision (a-h,o-z)
      dimension    q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)
      dimension  aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)
c
      return
      end
