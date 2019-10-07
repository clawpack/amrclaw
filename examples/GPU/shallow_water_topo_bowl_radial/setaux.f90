!     ============================================
subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
!     ============================================
!     
!     # set auxiliary arrays 
!     # variable coefficient acoustics
!     #  aux(1,i,j) =  topo in cell (i,j)

    implicit none
    integer, intent(in) :: mbc,mx,my,maux
    real(CLAW_REAL), intent(in) :: xlower,ylower,dx,dy
    real(CLAW_REAL), intent(inout) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer :: i,j
    real(CLAW_REAL) :: xcell, ycell
    real(CLAW_REAL), parameter :: zmin = 80.d0

    do i=1-mbc,mx+mbc
        xcell = xlower + (i-0.5d0)*dx
        do j=1-mbc,my+mbc
            ycell = ylower + (j-0.5d0)*dy
            aux(1,i,j) = -zmin + 0.01*(xcell**2 + ycell**2)
        enddo
    enddo

    return
end subroutine setaux
