
c
c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower,
     &                  dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
c
       dimension   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
       dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
c
c
      do 215 i=1,mx
        xcell = xlower + (i-0.5d0)*dx
        do 210 j=1,my
          ycell = ylower + (j-0.5d0)*dy
          do 205 k=1,mz
            zcell = zlower + (k-0.5d0)*dz
	          r = dsqrt(xcell**2 + ycell**2 + zcell**2)
	          rho = 1.d0 + 10.0d0*dexp(-20*(r-0.0d0)**2)
	          q(1,i,j,k) = rho
	          q(2,i,j,k) = 0.d0
	          q(3,i,j,k) = 0.d0
	          q(4,i,j,k) = 0.d0
	          q(5,i,j,k) = rho
205       continue
210     continue
215   continue

      return
      end
