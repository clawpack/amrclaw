! Set initial conditions for q.
subroutine qinit(meqn, mbc, mx, my, mz, xlower, ylower, zlower, dx, dy, dz, q, maux, aux)

    implicit none

    ! Input/Output
    integer, intent(in) :: meqn, mbc, mx, my, mz, maux
    real(kind=8), intent(in) :: xlower, ylower, zlower, dx, dy, dz
    real(kind=8), intent(inout) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)
    real(kind=8), intent(in) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)

    ! Local storage
    integer :: i, j, k
    real(kind=8) :: x, y, z, r, rho

    do i = 1, mx
        x = xlower + (i - 0.5d0) * dx
        do j = 1, my
            y = ylower + (j - 0.5d0) * dy
            do k = 1, mz
                z = zlower + (k - 0.5d0) * dz
                r = sqrt(x**2 + y**2 + z**2)
                rho = 1.d0 + 10.d0 * exp(-20.d0 * (r - 0.d0)**2)
                q(1, i, j, k) = rho
                q(2, i, j, k) = 0.d0
                q(3, i, j, k) = 0.d0
                q(4, i, j, k) = 0.d0
                q(5, i, j, k) = rho
            end do
        end do
    end do
end subroutine qinit
