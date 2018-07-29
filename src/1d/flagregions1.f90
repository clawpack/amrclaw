! ::::::::::::::::::::: flagregions ::::::::::::::::::::::::::::::::::
!
! Modify array of flagged points to respect minlevels and maxlevels
! specified by regions.
!
! Second version with outer loop on regions, should be faster.
!
! amrflags  = array to be flagged with either the value
!             DONTFLAG (no refinement needed)  or
!             DOFLAG   (refinement desired)    
!
! This routine may change flags only in cells that are (partially)
! covered by one or more regions. amrflags will be later modified
! by Richardson extrapolation and/or flag2refine routine, as requested,
! which will only add DOFLAG points to cells that are still UNSET
!
! If any part of a grid cell is covered by one or more regions, then
! refinement is *required* to at least the max of all region min_levels
! and is *allowed* to at most to the max of all region max_levels.
!
! Note that buffering is done *after* this, so additional cells may
! be refined in areas covered by regions that do not allow it!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine flagregions1(mx,mbuff,xlower,dx,level,t, &
                            amrflags)

    use regions_module
    use amr_module, only : DOFLAG, UNSET, DONTFLAG

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,level,mbuff
    real(kind=8), intent(in) :: xlower,dx,t
    
    ! Flagging
    real(kind=8),intent(inout) :: amrflags(1-mbuff:mx+mbuff)
    
    ! Locals
    integer :: i,m,i1,i2
    real(kind=8) :: x_low,x_hi, xupper
    integer, allocatable :: minlevel(:), maxlevel(:)

    allocate(minlevel(mx), maxlevel(mx))
    
    minlevel = 0
    maxlevel = 0

    xupper = xlower + mx*dx

    rloop: do m=1,num_regions
        if (t < regions(m)%t_low .or. t > regions(m)%t_hi) then
            cycle rloop  ! no intersection
        endif
        
        if (xlower >= regions(m)%x_hi .or. xupper <= regions(m)%x_low) then
            cycle rloop  ! no intersection
        else
            i1 = max(floor((regions(m)%x_low - xlower) / dx) + 1, 1)
            i2 = min(floor((regions(m)%x_hi -xlower) / dx) + 1, mx)
        endif

        do i=i1,i2
            minlevel(i) = max(minlevel(i), regions(m)%min_level)
            maxlevel(i) = max(maxlevel(i), regions(m)%max_level)
        enddo
    enddo rloop

    do i=1,mx
        if (minlevel(i) > maxlevel(i)) then
            write(6,*) '*** Error: this should never happen!'
            write(6,*) '*** minlevel > maxlevel in flagregions'
            stop
        endif

        if (maxlevel(i) /= 0) then
            ! this point lies in at least one region, so may need
            ! to modify the exisiting flag at this point...
            if (level < minlevel(i)) then
                ! Require refinement of this cell:
                amrflags(i) = DOFLAG
            else if (level >= maxlevel(i)) then
                ! Do not refine of this cell:
                amrflags(i) = DONTFLAG
            ! else leave amrflags(i) alone.
            endif
        endif

    enddo

end subroutine flagregions1
