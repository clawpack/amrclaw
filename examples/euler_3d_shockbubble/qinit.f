c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower,
     &                  dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
       dimension   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
       dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)


       double precision xc(1-mbc:mx+mbc)
       double precision xe(1-mbc:mx+mbc)

       double precision yc(1-mbc:my+mbc)
       double precision ye(1-mbc:my+mbc)

       double precision zc(1-mbc:mz+mbc)
       double precision ze(1-mbc:mz+mbc)

       common /cparam/ gamma
       common /comic/ qin(5),qout(5)
       common /cominf/ rinf,vinf,einf
       common /cdisc/ x0,y0,r0,radius

       do i = 1-mbc,mx+mbc
          xe(i) = xlower + (i - 1  )*dx
          xc(i) = xlower + (i - 0.5)*dx
       enddo

       do j = 1-mbc,my+mbc
          ye(j) = ylower + (j - 1.0)*dy
          yc(j) = ylower + (j - 0.5)*dy
       enddo

       do k = 1-mbc,mz+mbc
          ze(k) = zlower + (k - 1.0)*dz
          zc(k) = zlower + (k - 0.5)*dz
       enddo


       do k = 1,mz
          do j = 1,my
             do i = 1,mx
                if (zc(k) <= r0) then
                   radius = dsqrt(r0**2 - zc(k)**2)
                   call cellave(xe(i),ye(j),dx,dy,win)
                else
                   win = 0.0
                endif
                do m=1,meqn
                   q(m,i,j,k) = win*qin(m) + (1.d0-win)*qout(m)
                enddo

                if (xe(i) < 0.2) then
c                  # behind shock:
                   q(1,i,j,k) = rinf
                   q(2,i,j,k) = rinf*vinf
                   q(3,i,j,k) = 0.d0
                   q(4,i,j,k) = 0.d0
                   q(5,i,j,k) = einf
                endif
             enddo
          enddo
       enddo

       return
       end
