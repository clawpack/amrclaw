      subroutine setprob
      implicit double precision (a-h,o-z)
      common /param/  gamma,gamma1
c
c     # Euler equations

      gamma = 1.4d0
      gamma1 = gamma - 1.d0

      return
      end
