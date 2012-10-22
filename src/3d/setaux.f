c     ==================================================================
      subroutine setaux(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower,ylower,
     &                  zlower,dx,dy,dz,maux,aux)
c     ==================================================================
c
c     # set auxiliary arrays
c     # dummy routine when no auxiliary arrays
c
c
      implicit double precision (a-h,o-z)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux)
c
      do 10 k=1-mbc,mz+mbc
         do 10 j=1-mbc,my+mbc
            do 10 i=1-mbc,mx+mbc
               aux(i,j,k,1) = 1.d0
               aux(i,j,k,2) = 0.5d0
               aux(i,j,k,3) = 0.7d0
   10          continue
       return
       end
