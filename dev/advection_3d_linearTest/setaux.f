c     ==================================================================
      subroutine setaux(mbc,mx,my,mz,xlower,ylower,zlower,
     &                  dx,dy,dz,maux,aux)
c     ==================================================================
c
c     # set auxiliary arrays
c
c     # advection
c     #    aux(i,j,k,1) is u velocity on left face of cell
c     #    aux(i,j,k,2) is v velocity on bottom face of cell
c     #    aux(i,j,k,3) is w velocity on back face of cell
c
c
      implicit none

      integer  mx, my, mz, mbc, maux, i,j,k
      double precision xlower, ylower, zlower, dx, dy, dz
      double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      double precision ubar, vbar, wbar

      common /cparam/ ubar,vbar,wbar

c

      do  k = 1-mbc,mz+mbc
         do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
               aux(1,i,j,k) = ubar
               aux(2,i,j,k) = vbar
               aux(3,i,j,k) = wbar
            enddo
         enddo
      enddo

      return
      end
