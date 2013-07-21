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
      dimension  aux(maux,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &                                     1-mbc:maxmz+mbc)
c

       return
       end
