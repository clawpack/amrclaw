c
c
c
c     ==================================================================
      subroutine rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,maux,wave,s,amdq,apdq)
c     ==================================================================
c
c     # Riemann-solver for the 3d Burgers' equation
c     #    q_t  +  u*(.5*q^2)_x + v*(.5*q^2)_y + w*(.5*q^2)_z = 0
c     # where u,v,w are a given scalars, stored in the vector coeff
c     # that is set in setprob.f and passed in the common block comrp.
c
c       -----------------------------------------------------------
c
c     # solve Riemann problems along one slice of data.
c     # This data is along a slice in the x-direction if ixyz=1
c     #                               the y-direction if ixyz=2.
c     #                               the z-direction if ixyz=3.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, wave contains the waves, s the speeds,
c     # and amdq, apdq the left-going and right-going flux differences,
c     # respectively.  
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
      implicit real*8(a-h,o-z)
c
      dimension wave(meqn,mwaves,1-mbc:maxm+mbc)
      dimension    s(mwaves,1-mbc:maxm+mbc)
      dimension   ql(meqn,1-mbc:maxm+mbc)
      dimension   qr(meqn,1-mbc:maxm+mbc)
      dimension amdq(meqn,1-mbc:maxm+mbc)
      dimension apdq(meqn,1-mbc:maxm+mbc)
      dimension auxl(maux,1-mbc:maxm+mbc)
      dimension auxr(maux,1-mbc:maxm+mbc)
      logical, parameter :: efix = .true.
      common /comrp/ coeff(3)
c
c     # Set wave, speed, and flux differences:
c     ------------------------------------------
c
      do 30 i = 2-mbc, mx+mbc
         wave(1,1,i) = ql(1,i) - qr(1,i-1)
         s(1,i) = coeff(ixyz) * 0.5d0*(qr(1,i-1) + ql(1,i))

c        # The flux difference df = s*wave all goes in the downwind direction:
         amdq(1,i) = dmin1(s(1,i), 0.d0) * wave(1,1,i)
         apdq(1,i) = dmax1(s(1,i), 0.d0) * wave(1,1,i)

         if (efix) then
c           # entropy fix for transonic rarefactions:
            if (qr(1,i-1).lt.0.d0 .and. ql(1,i).gt.0.d0) then
               amdq(1,i) = - coeff(ixyz) * 0.5d0*qr(1,i-1)**2
               apdq(1,i) =   coeff(ixyz) * 0.5d0*ql(1,i)**2
               endif
            endif
   30    continue
c
      return
      end

