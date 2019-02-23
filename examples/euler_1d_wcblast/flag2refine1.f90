! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! -----------------------------------------------------------------------
! This version for wcblast problem uses only the jump in density to flag.
! -----------------------------------------------------------------------
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
! tolsp = tolerance specified by user in input file amr1ez.data, used in default
!         version of this routine as a tolerance for spatial differences.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine flag2refine1(mx,mbc,mbuff,meqn,maux,xlower,dx,t,level, &
                            tolsp,q,aux,amrflags)

    use regions_module
    use amr_module, only : DOFLAG, UNSET

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,dx,t,tolsp
    
    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc)
    
    ! Flagging
    real(kind=8),intent(inout) :: amrflags(1-mbuff:mx+mbuff)
    
    logical :: allowflag
    external allowflag

    ! Locals
    integer :: i,m
    real(kind=8) :: x_c,x_low,x_hi
    real(kind=8) :: dqi(meqn), dq(meqn)

    ! Don't initialize flags, since they were already 
    ! flagged by flagregions2
    ! amrflags = DONTFLAG
    
    ! Loop over interior points on this grid
    ! (i) grid cell is [x_low,x_hi], cell center at (x_c)
    ! This information is not needed for the default flagging based on
    ! undivided differences, but might be needed in a user's version.
    ! Note that if you want to refine only in certain space-time regions,
    ! it may be easiest to use the "regions" feature.  The flags set here or
    ! in the Richardson error estimator are potentially modified by the
    ! min_level and max_level specified in any regions.
        
    x_loop: do i = 1,mx
        !x_low = xlower + (i - 1) * dx
        !x_c = xlower + (i - 0.5d0) * dx
        !x_hi = xlower + i * dx

        ! -----------------------------------------------------------------
            ! Only check undivided differences if flag hasn't been set yet. 
            ! If flag == DONTFLAG then refinement is forbidden by a region, 
            ! if flag == DOFLAG checking is not needed
            if(amrflags(i) == UNSET) then
                dq = 0.d0
                dqi = abs(q(:,i+1) - q(:,i-1))
                dq = max(dq,dqi)

                ! check only undivided difference of density:
                do m=1,1
                    if (dq(m) > tolsp) then
                        amrflags(i) = DOFLAG
                        cycle x_loop
                    endif
                enddo
            endif

    enddo x_loop
   
end subroutine flag2refine1
