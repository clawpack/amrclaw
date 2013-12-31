subroutine setaux(mbc,mx,my,mz,xlower,ylower,zlower,dx,dy,dz,maux,aux,aux_copy_mask)
    ! set auxiliary arrays
    ! acoustics in a heterogeneous medium:
    ! aux(i,j,k,1) = impedance Z in (i,j) cell
    ! aux(i,j,k,2) = sound speed c in (i,j) cell

    implicit none

    integer, intent(in) :: mbc, mx, my, mz, maux
    real(kind=8), intent(in) :: xlower, ylower, zlower, dx, dy, dz
    real(kind=8), intent(in out) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    integer(kind=1), intent(in) :: aux_copy_mask(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    
    integer :: k, j, i
    real(kind=8) :: zcell, z1, c1, z2, c2

    common /comaux/ z1,c1,z2,c2


    do  k = 1-mbc,mz+mbc
        do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
                zcell = zlower + (k-0.5d0)*dz
                if (zcell .lt. 0.d0) then
                    aux(1,i,j,k) = z1
                    aux(2,i,j,k) = c1
                else
                    aux(1,i,j,k) = z2
                    aux(2,i,j,k) = c2
                endif
            enddo
        enddo
    enddo

end subroutine setaux