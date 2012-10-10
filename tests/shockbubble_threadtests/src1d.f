c
c
c =========================================================
      subroutine src1d(meqn,mbc,mx1d,q1d,maux,aux1d,t,dt)
c =========================================================
      implicit double precision (a-h,o-z)
      dimension   q1d(meqn,mx1d)
      dimension aux1d(maux,mx1d)
c
c
c     # This routine should be a simplified version of src2 
c     # which applies source terms for a 1-d slice of data along the
c     # edge of a grid.  This is called only from qad where the conservative
c     # fix-up is applied and is used to apply source terms over partial
c     # time steps to the coarse grid cell values used in solving Riemann 
c     # problems at the interface between coarse and fine grids.
c
c     # If the source terms depend only on q, it should be easy to 
c     # adapt src2 to create this routine, just loop over 1:mx1d.
c     # If the source terms are more complicated, it may not be easy.
c
c     # The code may work fine without applying source terms in this
c     # context, so using this dummy routine might be successful even when
c     # source terms are present. 
c
c
      dimension qstar(4)
      common /cparam/  gamma,gamma1
c
c     # 2-stage Runge-Kutta method
c
      dt2    = dt/2.d0
      press  = 0.d0
      ndim   = 2

      do 10 i=1,mx1d
         rad      = aux1d(1,i)
         rho      = q1d(1,i)
         u        = q1d(2,i)/q1d(1,i)
         v        = q1d(3,i)/q1d(1,i)
         press    = gamma1*(q1d(4,i) - 0.5d0*rho*(u**2 + v**2))

         if (rad.eq.0.d0) write(6,*) 'rad = 0 in src2'
         if (rho.eq.0.d0) write(6,*) 'rho = 0 in q'

         qstar(1) = q1d(1,i) - dt2*(ndim-1)/rad * q1d(3,i)
         qstar(2) = q1d(2,i) - dt2*(ndim-1)/rad * 
     &                          (rho*u*v)
         qstar(3) = q1d(3,i) - dt2*(ndim-1)/rad * 
     &                          (rho*v*v)
         qstar(4) = q1d(4,i) - dt2*(ndim-1)/rad * 
     &                          v*(q1d(4,i) + press)
c
c        # second stage
c
         rho      = qstar(1)
         u        = qstar(2)/qstar(1)
         v        = qstar(3)/qstar(1)
         press    = gamma1*(qstar(4) - 0.5d0*rho*(u**2 + v**2))
         if (rho.eq.0.d0) write(6,*) 'rho = 0 in qstar'

         q1d(1,i) = q1d(1,i) - dt*(ndim-1)/rad * qstar(3)
         q1d(2,i) = q1d(2,i) - dt*(ndim-1)/rad * 
     &                          (rho*u*v)
         q1d(3,i) = q1d(3,i) - dt*(ndim-1)/rad * 
     &                          (rho*v*v)
         q1d(4,i) = q1d(4,i) - dt*(ndim-1)/rad * 
     &                          v*(qstar(4) + press)
   10    continue

      return
      end
