c
c     =================================================
      function fdisc(x,y)
c     =================================================
      implicit double precision (a-h,o-z)
      common/cdisc/ x0,y0, r0, radius
c
c     # for computing cell averages for initial data that has a
c     # discontinuity along some curve.  fdisc should be negative to the
c     # left of the curve and positive to the right
c     # idisc specifies the nature of the discontinuity for two
c     # particular cases (a straight line and circle) but this routine
c     # can be modified for any other curve.
c

c     # circle of radius r0:


      fdisc = (x-x0)**2 + (y-y0)**2 - radius**2
c
      return
      end
