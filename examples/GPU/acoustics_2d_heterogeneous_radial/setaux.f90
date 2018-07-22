!     ============================================
subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
!     ============================================
!     
!     # set auxiliary arrays 
!     # variable coefficient acoustics
!     #  aux(1,i,j) =  density rho in cell (i,j)
!     #  aux(2,i,j) = bulk modulus in cell (i,j)
!     
!     # Piecewise constant medium with single interface at x=0
!     # Density and bulk modulus to left and right are set in problem_para_module.f90
    use problem_para_module, only: rho_l, rho_r, bulk_l, bulk_r

    implicit none
    integer, intent(in) :: mbc,mx,my,maux
    real(CLAW_REAL), intent(in) :: xlower,ylower,dx,dy
    real(CLAW_REAL), intent(inout) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer i,j
    double precision xcell

    do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
            xcell = xlower + (i-0.5d0)*dx
            if (xcell .lt. 0.0d0) then
                aux(1,i,j) = rho_l
                aux(2,i,j) = bulk_l
            else
                aux(1,i,j) = rho_r
                aux(2,i,j) = bulk_r
            endif
        enddo
    enddo

    return
end subroutine setaux
