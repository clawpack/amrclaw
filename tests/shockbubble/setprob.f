c     ==================
      subroutine setprob()
c     ==================

      implicit double precision (a-h,o-z)
      common /comic/ qin(5),qout(5)
      common /cparam/  gamma,gamma1
      common/cdisc/ x0,y0,alf,beta,r0,idisc
      common /cominf/ rinf,vinf,einf
c
c
C       open(unit=7,file='setprob.data',status='old',form='formatted')
	   call opendatafile(7,"setprob.data")
c      # set idisc for cellave routines (see function fdisc) 
       idisc = 2
c
       read(7,*) gamma
       gamma1 = gamma - 1.d0

c      # read center and radius of bubble:
       read(7,*) x0,y0,r0
c      # density in bubble:
       read(7,*) rhoin
c      # pressure behind shock:
       read(7,*) pinf 
c
c      # density outside bubble and pressure ahead of shock are fixed:
       rhoout = 1.d0
       pout   = 1.d0
       pin    = 1.d0

       qin(1) = rhoin
       qin(2) = 0.d0
       qin(3) = 0.d0
       qin(4) = pin/gamma1 
       qin(5) = 1.d0

       qout(1) = rhoout
       qout(2) = 0.d0
       qout(3) = 0.d0
       qout(4) = pout/gamma1 
       qout(5) = 0.d0
c
c     # Compute density and velocity behind shock from Hugoniot relations:
c     # ------------------------------------------------------------------

      rinf = ( gamma1 + (gamma+1)*pinf )/( (gamma+1) + gamma1*pinf )
      vinf = (1.0d0/sqrt(gamma)) * (pinf - 1.d0)/
     &       sqrt( 0.5*((gamma+1)/gamma) * pinf + 0.5*gamma1/gamma )
      einf = 0.5*rinf*vinf*vinf + pinf/gamma1

      write(6,601) pinf,rinf,vinf,einf
  601 format('pinf,rinf,vinf,einf:',/, 4e14.6)
      write(6,*) " ******** max1d = ",max1d
      return
      end
