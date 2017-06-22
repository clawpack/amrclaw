c     ============================================
      subroutine setaux(mbc,mx,xlower,dx,maux,aux)
c     ============================================
c     
c     # set auxiliary arrays 
c     # variable coefficient acoustics
c     #  aux(i,1) = impedance Z in i'th cell
c     #  aux(i,2) = sound speed c in i'th cell
c     
c     # Piecewise constant medium with single interface at x=0
c     # Density and sound speed to left and right are set in setprob.f
c

      implicit none

      integer, intent(in) :: mbc, mx, maux
      double precision, intent(in) :: xlower, dx
      double precision, intent(out) :: aux
      dimension aux(maux, 1-mbc:mx+mbc)

      common /comaux/ Zl, cl, Zr, cr
      double precision Zl, cl, Zr, cr

      integer i,ii
      double precision xcell

      do i=1-mbc,mx+mbc
         xcell = xlower + (i-0.5d0)*dx
         if (xcell .lt. 0.0d0) then
            aux(1,i) = Zl
            aux(2,i) = cl
         else
            aux(1,i) = Zr
            aux(2,i) = cr
         endif
      enddo

      return
      end
