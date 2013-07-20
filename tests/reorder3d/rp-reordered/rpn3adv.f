c
c
c
c     ==================================================================
      subroutine rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,maux,wave,s,amdq,apdq)
c     ==================================================================
c
c     # Riemann-solver for the advection equation
c     #    q_t  +  u*q_x + v*q_y + w*q_z = 0
c     # where u and v are a given velocity field.
c
c       -----------------------------------------------------------
c     # In advective form, with interface velocities specified in
c     # the auxiliary variable
c     # auxl(i,ma) contains auxiliary data for cells along this slice,
c     #   ma=1:   u-velocity at left edge of cell 
c     #   ma=2:   v-velocity at bottom edge of cell
c     #   ma=3:   w-velocity at back edge of cell 
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
c     # respectively.  Note that in this advective form, the sum of
c     # amdq and apdq is not equal to a difference of fluxes except in the
c     # case of constant velocities.
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
      implicit real*8(a-h,o-z)
c
      dimension wave(meqn,mwaves,1-mbc:maxm+mbc)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(meqn,1-mbc:maxm+mbc)
      dimension   qr(meqn,1-mbc:maxm+mbc)
      dimension amdq(meqn,1-mbc:maxm+mbc)
      dimension apdq(meqn,1-mbc:maxm+mbc)
      dimension auxl(maux,1-mbc:maxm+mbc)
      dimension auxr(maux,1-mbc:maxm+mbc)
c
c
c     # Set wave, speed, and flux differences:
c     ------------------------------------------
c
      do 30 i = 2-mbc, mx+mbc
         wave(1,1,i) = ql(1,i) - qr(1,i-1)
         s(i,1) = auxl(ixyz,i)
c        # The flux difference df = s*wave all goes in the downwind direction:
         amdq(1,i) = dmin1(auxl(ixyz,i), 0.d0) * wave(1,1,i)
         apdq(1,i) = dmax1(auxl(ixyz,i), 0.d0) * wave(1,1,i)
   30    continue
c
      return
      end

