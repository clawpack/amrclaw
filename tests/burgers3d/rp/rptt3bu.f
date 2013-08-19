c
c
c     ==================================================================
      subroutine rptt3(ixyz,icoor,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,maux,imp,impt,bsasdq,
     &                  cmbsasdq,cpbsasdq)
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
c
c     #    bsasdq is an array of flux differences that result from a
c     #         transverse splitting (a previous call to rpt3).  
c     #         This stands for B^* A^* \Dq but could represent any of 
c     #         6 possibilities, e.g.  C^* B^* \Dq, as specified by ixyz
c     #         and icoor (see below).
c     #         Moreover, each * represents either + or -, as specified by
c     #         imp and impt.
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
c     #        ixyz=1, icoor=3 means bsasdq=B^*A^*\Dq, and should be 
c     #                        split in z into 
c     #                           cmbsasdq = C^-B^*A^*\Dq,
c     #                           cpbsasdq = C^+B^*A^*\Dq.
c     #
c     #        ixyz=2, icoor=3 means bsasdq=C^*B^*\Dq, and should be
c     #                        split in x into 
c     #                           cmbsasdq = A^-C^*B^*\Dq,
c     #                           cpbsasdq = A^+C^*B^*\Dq.
c
c     #    The parameters imp and impt are generally needed only if aux
c     #    arrays are being used, in order to access the appropriate
c     #    variable coefficients.

c
      implicit real*8(a-h,o-z)
      dimension       ql(meqn,1-mbc:maxm+mbc)
      dimension       qr(meqn,1-mbc:maxm+mbc)
      dimension   bsasdq(meqn,1-mbc:maxm+mbc)
      dimension cmbsasdq(meqn,1-mbc:maxm+mbc)
      dimension cpbsasdq(meqn,1-mbc:maxm+mbc)
      dimension     aux1(maux,1-mbc:maxm+mbc,3)
      dimension     aux2(maux,1-mbc:maxm+mbc,3)
      dimension     aux3(maux,1-mbc:maxm+mbc,3)
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
          sc = coeff(iuvw) * 0.5d0*(qr(1,i-1) + ql(1,i))
          cmbsasdq(1,i) = dmin1(sc, 0.d0) * bsasdq(1,i)
          cpbsasdq(1,i) = dmax1(sc, 0.d0) * bsasdq(1,i)
   10    continue
c
      return
      end

