subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux,aux_copy_mask)
    ! Set auxiliary arrays 
    ! aux(1,i,j) is edge velocity at "left" boundary of grid point (i,j)
    ! aux(2,i,j) is edge velocity at "bottom" boundary of grid point (i,j)
    ! aux(3,i,j) is kappa if a mapped grid is used.
 
    implicit none
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(in out) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    integer(kind=1), intent(in) :: aux_copy_mask(1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Local storage
    integer :: i, j
    real(kind=8) :: xll, yll

    ! Stream function interface
    interface
        real(kind=8) function psi(x,y)
            implicit none
            real(kind=8), intent(in) :: x, y
        end function psi
    end interface

    ! constant velocities which are used if tperiod=0 is specified
    ! in setprob.data

    do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc

            ! coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

            ! difference stream function psi to get normal velocities:
            aux(1,i,j) = -(psi(xll, yll+dy) - psi(xll,yll)) / dy
            aux(2,i,j) =  (psi(xll+dx, yll) - psi(xll,yll)) / dx
        end do
    end do

end subroutine setaux
