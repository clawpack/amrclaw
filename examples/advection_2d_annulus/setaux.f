c     ============================================
      subroutine setaux(mbc,mx,my,xlow,ylow,dxc,dyc,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays for advection on a curvilinear grid
c
c     # on input, (xc(i),yc(j)) gives uniformly spaced computational grid.
c     # on output, 
c     #   aux(1,i,j) is edge velocity at "left" boundary of grid point (i,j)
c     #   aux(2,i,j) is edge velocity at "bottom" boundary of grid point (i,j)
c     #   aux(3,i,j) = kappa  is ratio of cell area to dxc*dyc
c     
      use amr_module, only: NEEDS_TO_BE_SET, xlower, ylower
      implicit double precision (a-h,o-z)
      double precision, intent(inout), 
     &              dimension(3, 1-mbc:mx+mbc,1-mbc:my+mbc) :: aux
      dimension xccorn(5),yccorn(5),xpcorn(5),ypcorn(5)
c
          ilo = floor((xlow - xlower + .05*dxc)/dxc)
          jlo = floor((ylow - ylower + .05*dyc)/dyc)
c
      do 20 j=1-mbc,my+mbc
         do 20 i=1-mbc,mx+mbc
c
            if (aux(1,i,j) .ne. NEEDS_TO_BE_SET) cycle
c
c           # computational points (xc,yc) are mapped to physical
c           # coordinates (xp,yp) by mapc2p:
c
c           # lower left corner:
!           xccorn(1) = xlower + (i-1)*dxc orig version, with xlower !           passed in
!           yccorn(1) = ylower + (j-1)*dyc
!           xccorn(1) = xlow + (i-1)*dxc   ! changed names of args, so !           this is equiv.
!           yccorn(1) = ylow + (j-1)*dyc   ! next should be w/o roundoff
            xccorn(1) = xlower + (ilo + i-1)*dxc
            yccorn(1) = ylower + (jlo + j-1)*dyc
            call mapc2p(xccorn(1),yccorn(1),xpcorn(1),ypcorn(1))

c           # upper left corner:
            xccorn(2) = xccorn(1)
            yccorn(2) = yccorn(1) + dyc
            call mapc2p(xccorn(2),yccorn(2),xpcorn(2),ypcorn(2))
c
c           # upper right corner:
            xccorn(3) = xccorn(1) + dxc
            yccorn(3) = yccorn(1) + dyc
            call mapc2p(xccorn(3),yccorn(3),xpcorn(3),ypcorn(3))
c
c           # lower right corner:
            xccorn(4) = xccorn(1) + dxc
            yccorn(4) = yccorn(1)
            call mapc2p(xccorn(4),yccorn(4),xpcorn(4),ypcorn(4))
c
c
c           # compute edge velocities by differencing stream function:
c
            aux(1,i,j) = (stream(xpcorn(2),ypcorn(2))
     &                    - stream(xpcorn(1),ypcorn(1)))/ dyc
c
            aux(2,i,j) = -(stream(xpcorn(4),ypcorn(4))
     &                    - stream(xpcorn(1),ypcorn(1)))/ dxc

c
c
c           # compute area of physical cell from four corners:

            xpcorn(5) = xpcorn(1)
            ypcorn(5) = ypcorn(1)
            area = 0.d0
            do ic=1,4
              area = area + 0.5d0 * (ypcorn(ic)+ypcorn(ic+1)) *
     &               (xpcorn(ic+1)-xpcorn(ic))
            enddo
            
            aux(3,i,j) = area / (dxc*dyc)
c
   20  continue
c
       return

       end
