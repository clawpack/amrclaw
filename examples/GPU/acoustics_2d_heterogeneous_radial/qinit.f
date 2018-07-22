
c
c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c     # Acoustics with smooth radially symmetric profile to test accuracy
c
       implicit none
    
       integer :: meqn, mbc, mx, my, maux
       real(CLAW_REAL) :: xlower, ylower, dx, dy
       real(CLAW_REAL) :: q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
       real(CLAW_REAL) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
       real(CLAW_REAL) :: x0, y0, xcell, ycell
       real(CLAW_REAL) :: r, width, rad, pressure, pi
       integer :: i,j
c
       pi = 4.d0*atan(1.d0)
       width = 0.1d0
       rad = 0.25d0
       x0 = -0.5
       y0 = 0.0

       do 20 i=1,mx
          xcell = xlower + (i-0.5d0)*dx
          do 20 j=1,my
             ycell = ylower + (j-0.5d0)*dy
             r = sqrt((xcell-x0)**2 + (ycell-y0)**2)

             if (abs(r-rad) .le. width) then
                 pressure = 1.d0 + cos(pi*(r - rad)/width)
             else
                 pressure = 0.d0
             endif
             q(1,i,j) = pressure
             q(2,i,j) = 0.d0
             q(3,i,j) = 0.d0
  20         continue

       return
       end
