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
            aux(1,i,j) = -1.d0
        enddo
    enddo

    return
end subroutine setaux
