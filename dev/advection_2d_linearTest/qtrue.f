

c     ==================================
      double precision function qtrue(x,y,t)
c     ==================================
      implicit double precision (a-h,o-z)
      common /cparam/ ubar,vbar
      
      x0 = x - ubar*t
      y0 = y - vbar*t

c     # evaluate desired initial data at (x0,y0):

!--      r = dsqrt((x0+0.2d0)**2 + (y0-0.4d0)**2)
!--      if (r <= 0.3d0) then
!--          qtrue = 1.d0
!--        else
!--          qtrue = 0.d0
!--        endif
!--
!--      r = dsqrt((x0-0.3d0)**2 + (y0-0.1d0)**2)
!--      qtrue = qtrue - dexp(-15.d0*r**2)


c linear test function to see if ghost cell interpolation exact
      qtrue = 2.d0*x0 + 1.0d0*y0
      !qtrue = 1.d0*x0 + 1.0d0*y0

      return
      end
