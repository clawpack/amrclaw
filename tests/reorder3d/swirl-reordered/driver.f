      program driver
c
c  Generic driver routine for claw3
c
c  Author: Donna Calhoun
c  Date : 3/22/02
c  Version : To be used with Clawpack 4.0
c
      implicit double precision (a-h,o-z)

c     # set parameters for maximum array sizes used in declarations
c     # these must be increased for larger problems.
c
c
      parameter (maxmx =   100)
      parameter (maxmy =   100)
      parameter (maxmz =   100)
      parameter (mwork =  2132560)

      parameter (mbc = 2)
      parameter (meqn = 1)
      parameter (mwaves = 1)
      parameter (maux = 3)

      dimension
     & q(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)

      dimension
     & aux(maux,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)

      dimension mthlim(mwaves)
      dimension work(mwork)
c
      call claw3ez(maxmx,maxmy,maxmz,meqn,mwaves,mbc,maux,mwork,mthlim,
     &           q,work,aux)

      stop
      end

