
c
c
c
c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
       implicit double precision (a-h,o-z)
       dimension q(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
       common /comic/ qin(5),qout(5)
       common /cominf/ rinf,vinf,einf
c
c
       do 50 i=1,mx
          xlow = xlower + (i-1)*dx
          do 20 j=1,my
c            # set (xlow,ylow) to lower left corner of grid cell:
             ylow = ylower + (j-1)*dy
             call cellave(xlow,ylow,dx,dy,win)
c            # win is now the fraction of the cell that lies inside the circle
             do 10 m=1,meqn
                q(m,i,j) = win*qin(m) + (1.d0-win)*qout(m)
  10            continue
  20         continue
         if (xlow .lt. 0.2d0) then
c           # behind shock:
            do 30 j=1,my
               q(1,i,j) = rinf
               q(2,i,j) = rinf*vinf
               q(3,i,j) = 0.d0
               q(4,i,j) = einf
               q(5,i,j) = 0.d0
   30          continue
            end if
c        
c        if (xlow .lt. 0.5d0) then
c           # to give two different values of tracer in bubble
c           # for better visualization of motion:
c           do 40 j=1,my
c              q(5,i,j) = 2.d0*q(5,i,j)
c  40          continue
c           end if

   50    continue
       return
       end
