! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! Modified flag2refine file to use adjoint-flagging.
!
! Default version computes spatial difference dq in each direction and
! for each component of q and flags any point where this is greater than
! the tolerance tolsp.  This is consistent with what the routine errsp did in
! earlier versions of amrclaw (4.2 and before).
!
! This routine can be copied to an application directory and modified to
! implement some other desired refinement criterion.
!
! Points may also be flagged for refining based on a Richardson estimate
! of the error, obtained by comparing solutions on the current grid and a
! coarsened grid.  Points are flagged if the estimated error is larger than
! the parameter tol in amr2ez.data, provided flag_richardson is .true.,
! otherwise the coarsening and Richardson estimation is not performed!  
! Points are flagged via Richardson in a separate routine.
!
! Once points are flagged via this routine and/or Richardson, the subroutine
! flagregions is applied to check each point against the min_level and
! max_level of refinement specified in any "region" set by the user.
! So flags set here might be over-ruled by region constraints.
!
!    q   = grid values including ghost cells (bndry vals at specified
!          time have already been set, so can use ghost cell values too)
!
!  aux   = aux array on this grid patch
!
! amrflags  = array to be flagged with either the value
!             DONTFLAG (no refinement needed)  or
!             DOFLAG   (refinement desired)    
!
! tolsp = tolerance specified by user in input file amr2ez.data, used in default
!         version of this routine as a tolerance for spatial differences.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine flag2refine2(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,dx,dy,t,level, &
                            tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)

    use innerprod_module, only: calculate_innerproduct
    use adjoint_module, only: totnum_adjoints, innerprod_index
    use adjoint_module, only: adjoints, trange_start, trange_final

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp
    
    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: aux_temp(maux,1:mx,1:my)

    ! Flagging
    real(kind=8),intent(inout) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)
    real(kind=8), intent(in) :: DONTFLAG
    real(kind=8), intent(in) :: DOFLAG
    
    logical :: allowflag
    external allowflag

    ! Locals
    integer :: i, j, r
    logical :: mask_selecta(totnum_adjoints)

    ! Initialize flags
    amrflags = DONTFLAG
    aux(innerprod_index,:,:) = 0.0
    mask_selecta = .false.

    ! Loop over adjoint snapshots
    do r=1,totnum_adjoints
        if ((t+adjoints(r)%time) >= trange_start .and. &
          (t+adjoints(r)%time) <= trange_final) then
            mask_selecta(r) = .true.
        endif
    enddo

    do r=1,totnum_adjoints-1
        if((.not. mask_selecta(r)) .and. &
          (mask_selecta(r+1))) then
            mask_selecta(r) = .true.
            exit
        endif
    enddo

    do r=totnum_adjoints,2,-1
        if((.not. mask_selecta(r)) .and. &
          (mask_selecta(r-1))) then
            mask_selecta(r) = .true.
            exit
        endif
    enddo

    ! Loop over adjoint snapshots
    aloop: do r=1,totnum_adjoints

        ! Consider only snapshots that are within the desired time range
        if (mask_selecta(r)) then
            ! Calculate inner product with current snapshot
            aux_temp(innerprod_index,:,:) = &
                calculate_innerproduct(t,q,r,mx,my,xlower,ylower,dx,dy,meqn,mbc)

            ! Save max inner product
            do i=1,mx
                do j = 1,my
                    aux(innerprod_index,i,j) = &
                        max(aux(innerprod_index,i,j), &
                        aux_temp(innerprod_index,i,j))
                enddo
            enddo

        endif

    enddo aloop

    ! Flag locations that need refining
    y_loop: do j = 1,my
      x_loop: do i = 1,mx
        if (aux(innerprod_index,i,j) > tolsp) then
            amrflags(i,j) = DOFLAG
            cycle x_loop
        endif
      enddo x_loop
    enddo y_loop

end subroutine flag2refine2
