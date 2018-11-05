!     =====================================================
subroutine qinit(meqn,mbc,mx,my,xlower,ylower, &
        dx,dy,q,maux,aux)
    !     =====================================================
    !
    !     # Set initial conditions for q.
    implicit none

    integer :: meqn, mbc, mx, my, maux
    real(CLAW_REAL) :: xlower, ylower, dx, dy
    real(CLAW_REAL) :: q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(CLAW_REAL) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(CLAW_REAL) :: xcell, ycell
    real(CLAW_REAL) :: topo, ze, eta, depth
    integer :: i,j

    real(CLAW_REAL), parameter :: zmin = 80.d0

    do i=1,mx
        xcell = xlower + (i-0.5d0)*dx
        do j=1,my
            ycell = ylower + (j-0.5d0)*dy

            topo = -zmin + 0.01*(xcell**2 + ycell**2)
            ze = -(xcell**2 + ycell**2)/10.d0
            if (ze > -10.d0) then
                eta = 40.d0*exp(ze)
            else
                eta = 0.d0
            endif
            depth = eta-topo
            if (depth <0.d0) then
                depth = 0.d0
            endif
            q(1,i,j) = depth
            q(2,i,j) = 0.d0
            q(3,i,j) = 0.d0
        enddo
    enddo
    return
end
