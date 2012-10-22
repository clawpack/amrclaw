c     ============================================
      subroutine b4step3(maxmx,maxmy,maxmz,mbc,mx,my,mz,meqn,q,
     &            xlower,ylower,zlower,dx,dy,dz,t,dt,maux,aux)
c     ============================================
c
c     # called from claw3 before each call to step3.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.

c
c     # dummy routine 
c
c     
      implicit double precision (a-h,o-z)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 
     &               1-mbc:maxmz+mbc, maux)
c
      return
      end
