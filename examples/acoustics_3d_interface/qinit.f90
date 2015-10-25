! =====================================================
!  Set initial conditions for q.
! =====================================================
subroutine qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,q,maux,aux)

    implicit none

    ! Input
    integer, intent(in) :: meqn, mbc, mx, my, mz, maux
    real(kind=8), intent(in) :: xlower, ylower, zlower, dx, dy, dz
    real(kind=8), intent(in out) :: q(meqn, 1-mbc:mx+mbc,    &
                                            1-mbc:my+mbc,    &
                                            1-mbc:mz+mbc)
    real(kind=8), intent(in) :: aux(maux, 1-mbc:mx+mbc,      &
                                          1-mbc:my+mbc,      &
                                          1-mbc:mz+mbc)

    ! Local storage
    integer :: k, j, i
    real(kind=8) :: zcell

    q = 0.d0

    ! Jump discontinuity aligned with jump in z,c:
    do k = 1 - mbc , mz + mbc
        do j = 1 - mbc, my + mbc
            do i = 1 - mbc, mx + mbc
                zcell = zlower + (k - 0.5d0) * dz
                if (zcell < 0.d0) then
                    q(1,i,j,k) = 1.d0
                else
                    q(1,i,j,k) = 0.d0
                endif
            enddo
        enddo
    enddo

end subroutine qinit