

c     ==================================
      double precision function qtrue(x,y,z,t)
c     ==================================
      implicit double precision (a-h,o-z)
      common /cparam/ ubar,vbar,wbar
      
      x0 = x - ubar*t
      y0 = y - vbar*t
      z0 = z - wbar*t

c     # evaluate desired initial data at (x0,y0,z0):


c linear test function to see if ghost cell interpolation exact
      qtrue = 1.d0*x0 + 2.0d0*y0 + 3.d0*z0
      !qtrue = 1.d0*x0 
      !if (x0 .lt. xupper/2) then
      !    qtrue = 1.d0
      !else
      !    qtrue = 2.d0
      !endif

      return
      end
