c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
       dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)

       do 20 i=1,mx
          xcell = xlower + (i-0.5d0)*dx

          do 20 j=1,my
             ycell = ylower + (j-0.5d0)*dy

             if (dabs(xcell-3.5d0) .le. dx
     &            .and. dabs(ycell-0.5d0) .le. dy) then
                 pressure = 2.d0
             else
                 pressure = 0.d0
             endif
             q(1,i,j) = pressure
             q(2,i,j) = 0.d0
             q(3,i,j) = 0.d0
  20         continue
       return
       end
