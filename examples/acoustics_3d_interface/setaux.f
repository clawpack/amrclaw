c     ==================================================================
      subroutine setaux(mbc,mx,my,mz,xlower,ylower,zlower,
     &                  dx,dy,dz,maux,aux)
c     ==================================================================
c
c     # set auxiliary arrays
c
c     # acoustics in a heterogeneous medium:
c     #  aux(i,j,k,1) = impedance Z in (i,j) cell
c     #  aux(i,j,k,2) = sound speed c in (i,j) cell
c
c
      implicit double precision (a-h,o-z)
      double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      common /comaux/ z1,c1,z2,c2


      do  k = 1-mbc,mz+mbc
        do j = 1-mbc,my+mbc
           do i = 1-mbc,mx+mbc
               zcell = zlower + (k-0.5d0)*dz
               if (zcell .lt. 0.d0) then
                  aux(1,i,j,k) = z1
                  aux(2,i,j,k) = c1
                else
                  aux(1,i,j,k) = z2
                  aux(2,i,j,k) = c2
                endif
           enddo
        enddo
      enddo

      return
      end

