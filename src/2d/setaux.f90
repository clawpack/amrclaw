subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)

    ! Called at start of computation before calling qinit, and
    ! when AMR is used, also called every time a new grid patch is created.
    ! Use to set auxiliary arrays aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc).
    ! Note that ghost cell values may need to be set if the aux arrays
    ! are used by the Riemann solver(s).
    !
    ! This default version does nothing, unless using adjoint-flagging.

    use amr_module, only : NEEDS_TO_BE_SET
    use adjoint_module, only : adjoint_flagging,innerprod_index
 
    implicit none
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer :: jj, ii

    if(adjoint_flagging) then
        ! If a new grid has been created, but hadn't been flagged
        ! set innerproduct to zero.
        do jj=1-mbc,my+mbc
            do ii=1-mbc,mx+mbc
                if (aux(1,ii,jj) .eq. NEEDS_TO_BE_SET) then
                    aux(innerprod_index,ii,jj) = 0.d0
                endif
            enddo
        enddo
    endif

end subroutine setaux
