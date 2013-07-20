c
c
c     ==================================================================
      subroutine rpt3(ixyz,icoor,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,maux,imp,asdq,
     &                  bmasdq,bpasdq)
c     ==================================================================
c
c     # Riemann solver in the transverse direction for the
c     # advection equations.
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
c     #    variable coefficients:

c     #    imp = 1 if asdq = A^- \Dq,  the left-going flux difference
c     #          2 if asdq = A^+ \Dq, the right-going flux difference
c
c     #    aux2(:,:,2) is a 1d slice of the aux array along the row
c     #                 where the data ql, qr lie.   
c     #    aux1(:,:,2) and aux3(:,:,2) are neighboring rows in the 
c     #                 y-like direction
c     #    aux2(:,:,1) and aux2(:,:,3) are neighboring rows in the 
c     #                z-like direction
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
c
c
c     # set iuvw = 1 for u, 2 for v, 3 for w component of velocity
c     # depending on transverse direction:
      iuvw = ixyz + icoor - 1
      if (iuvw.gt.3) iuvw = iuvw-3
c
      do 10 i=2-mbc,mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
         if (icoor == 2) then  !! transverse dir. is y-like direction
            bmasdq(1,i) = dmin1(aux2(iuvw,i1,2), 0.d0)*asdq(1,i)
            bpasdq(1,i) = dmax1(aux3(iuvw,i1,2), 0.d0)*asdq(1,i)
         else !! icoor == 3  !! transverse dir. is z-like direction
            !! quanities split into cmasdq and cpasdq
            bmasdq(1,i) = dmin1(aux2(iuvw,i1,2),0.d0)*asdq(1,i)
            bpasdq(1,i) = dmax1(aux2(iuvw,i1,3),0.d0)*asdq(1,i)
         endif
   10    continue
c
      return
      end

