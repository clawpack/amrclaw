
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)
c     =====================================================
c
c     # Solve Riemann problems for the 2D hyperbolic problem.
c
c     # This is a wrapper for a 1d Riemann solver rp1 so that it
c     # can be used with AMR by using the 2d amrclaw routines with my=1.
c
      implicit double precision (a-h,o-z)
c
      dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension  apdq(meqn, 1-mbc:maxm+mbc)
      dimension  amdq(meqn, 1-mbc:maxm+mbc)
c
      if (ixy.eq.2) then
c        write(6,*) '*** Error, this Riemann solver should only be used'
c        write(6,*) '*** for a 1d problem with my=1'
c        write(6,*) '*** ixy=2 and my = ',mx
c        stop
         endif

      call rp1(maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)

      return
      end
