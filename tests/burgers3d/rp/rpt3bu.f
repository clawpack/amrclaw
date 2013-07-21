c
c
c     ==================================================================
      subroutine rpt3(ixyz,icoor,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,maux,imp,asdq,
     &                  bmasdq,bpasdq)
c     ==================================================================
c
c     # Riemann solver in the transverse direction for the 3d Burgers' equation
c     #    q_t  +  u*(.5*q^2)_x + v*(.5*q^2)_y + w*(.5*q^2)_z = 0
c     # where u,v,w are a given scalars, stored in the vector coeff
c     # that is set in setprob.f and passed in the common block comrp.
c     #
c     # On input,
c
c     #    ql,qr is the data along some one-dimensional slice, as in rpn3
c     #         This slice is
c     #             in the x-direction if ixyz=1,
c     #             in the y-direction if ixyz=2, or
c     #             in the z-direction if ixyz=3.
c     #    asdq is an array of flux differences (A^*\Dq).
c     #         asdq(i,:) is the flux difference propagating away from
c     #         the interface between cells i-1 and i.
c     #    Note that asdq represents B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.
c
c     #    ixyz indicates the direction of the original Riemann solve,
c     #         called the x-like direction in the table below:
c
c     #               x-like direction   y-like direction   z-like direction
c     #      ixyz=1:        x                  y                  z         
c     #      ixyz=2:        y                  z                  x         
c     #      ixyz=3:        z                  x                  y         
c
c     #    icoor indicates direction in which the transverse solve should 
c     #         be performed.
c     #      icoor=2: split in the y-like direction.
c     #      icoor=3: split in the z-like direction.
c
c     #    For example,
c     #      ixyz=1, icoor=2 means asdq=A^*\Dq, and should be split in y into
c     #                        bmasdq = B^-A^*\Dq,
c     #                        bpasdq = B^+A^*\Dq.
c     #
c     #      ixyz=2, icoor=2 means asdq=B^*\Dq, and should be split in z into 
c     #                        bmasdq = C^-B^*\Dq,
c     #                        bpasdq = C^+B^*\Dq.
c
c     #    The parameter imp is generally needed only if aux
c     #    arrays are being used, in order to access the appropriate
c     #    variable coefficients.
c
c
c
      implicit real*8(a-h,o-z)
      dimension     ql(meqn,1-mbc:maxm+mbc)
      dimension     qr(meqn,1-mbc:maxm+mbc)
      dimension   asdq(meqn,1-mbc:maxm+mbc)
      dimension bmasdq(meqn,1-mbc:maxm+mbc)
      dimension bpasdq(meqn,1-mbc:maxm+mbc)
      dimension   aux1(maux,1-mbc:maxm+mbc,3)
      dimension   aux2(maux,1-mbc:maxm+mbc,3)
      dimension   aux3(maux,1-mbc:maxm+mbc,3)
      common /comrp/ coeff(3)
c
c
c     # set iuvw = 1 for u, 2 for v, 3 for w component of velocity
c     # depending on transverse direction:
      iuvw = ixyz + icoor - 1
      if (iuvw.gt.3) iuvw = iuvw-3
c
c     # transverse wave goes up or down with speed based on data for
c     # original normal Riemann problem.  

      do 10 i=2-mbc,mx+mbc
          sb = coeff(iuvw) * 0.5d0*(qr(1,i-1) + ql(1,i))
          bmasdq(1,i) = dmin1(sb, 0.d0) * asdq(1,i)
          bpasdq(1,i) = dmax1(sb, 0.d0) * asdq(1,i)
   10    continue
c
      return
      end

