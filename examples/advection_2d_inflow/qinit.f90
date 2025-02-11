subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    ! Set initial conditions for q.
    ! Sample scalar equation with data that is piecewise constant with
    ! q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
    !     0.1  otherwise

    integer, intent(in) :: meqn, mbc, mx, my, maux
    real(kind=8), intent(in) :: xlower, ylower, dx, dy
    real(kind=8), intent(in out) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in out) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)

    real(kind=8), external :: qtrue

    ! set concentration profile
    do j=1,my
        yj = ylower + (j-0.5d0)*dy
        do i=1,mx
            xi = xlower + (i-0.5d0)*dx
            q(1,i,j) = qtrue(xi,yj,0.d0)
        end do
    end do

end subroutine qinit