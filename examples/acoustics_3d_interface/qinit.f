
c
c
c
c     =====================================================
      subroutine qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower,
     &                 dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
c
      dimension   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
c
c     # set initial data
c     ---------------------------
c
      q = 0.d0

c     # jump discontinuity aligned with jump in z,c:

      do  k = 1-mbc,mz+mbc
	    do j = 1-mbc,my+mbc
	       do i = 1-mbc,mx+mbc
               zcell = zlower + (k-0.5d0)*dz
               if (zcell .lt. 0.d0) then
                  q(1,i,j,k) = 1.d0
                else
                  q(1,i,j,k) = 0.d0
                endif
           enddo
        enddo
      enddo


      return
      end
