
c
c
c
c     =====================================================
       subroutine qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &                   xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
c
       dimension
     & q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)
       dimension x(1-mbc:maxmx+mbc)
       dimension 
     & aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)
c
c     # set concentration profile
c     ---------------------------
c

       do i = 1-mbc,mx+mbc
          x(i) = xlower + (i-0.5d0)*dx
       enddo

       do i = 1,mx
          if (x(i) .lt.  0.5d0) then
             do j = 1,my
                do k = 1,mz
                   q(1,i,j,k) = 1.d0
                enddo
             enddo
          else
             do j = 1,my
                do k = 1,mz
                   q(1,i,j,k) = 0.d0
                enddo
            enddo
         endif
      enddo

      return
      end

