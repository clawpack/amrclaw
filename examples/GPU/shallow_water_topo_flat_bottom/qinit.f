
c
c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
       implicit none
    
       integer :: meqn, mbc, mx, my, maux
       real(CLAW_REAL) :: xlower, ylower, dx, dy
       real(CLAW_REAL) :: q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
       real(CLAW_REAL) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
       real(CLAW_REAL) :: x0, y0, xcell, ycell
       real(CLAW_REAL) :: r, rad, depth, radius
       integer :: i,j

       radius = 0.5d0
       x0 = 0.d0
       y0 = 0.d0

       do 20 i=1,mx
          xcell = xlower + (i-0.5d0)*dx
          do 20 j=1,my
             ycell = ylower + (j-0.5d0)*dy
             r = sqrt((xcell-x0)**2 + (ycell-y0)**2)

             if (r .le. radius) then
                 depth = 2.d0
             else
                 depth = 1.d0
             endif
             q(1,i,j) = depth
             q(2,i,j) = 0.d0
             q(3,i,j) = 0.d0
  20         continue

       return
       end
