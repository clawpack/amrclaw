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

    integer i,j
    double precision xcell

    do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
!             xcell = xlower + (i-0.5d0)*dx
!             if (xcell .lt. 0.0d0) then
!                 aux(1,i,j) = rho_l
!                 aux(2,i,j) = bulk_l
!             else
!                 aux(1,i,j) = rho_r
!                 aux(2,i,j) = bulk_r
!             endif
            aux(1,i,j) = -1.d0
        enddo
    enddo

    return
end subroutine setaux
