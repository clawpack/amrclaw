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
! On input, amrflags is already set by flagger routine using 
! Richardson extrapolation and/or flag2refine routine, as requested.
! This routine may change flags only in cells that are (partially)
! covered by one or more regions.
!
! If any part of a grid cell is covered by one or more regions, then
! refinement is *required* to at least the max of all region min_levels
! and is *allowed* to at most to the max of all region max_levels.
!
! Note that buffering is done *after* this, so additional cells may
! be refined in areas covered by regions that do not allow it!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine flagregions2(mx,my,mbuff,xlower,ylower,dx,dy,level,t, &
                            amrflags,DONTFLAG,DOFLAG)

    use regions_module

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,level,mbuff
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t
    
    ! Flagging
    real(kind=8),intent(inout) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)
    real(kind=8), intent(in) :: DONTFLAG
    real(kind=8), intent(in) :: DOFLAG
    
    ! Locals
    integer :: i,j,m,i1,i2,j1,j2
    real(kind=8) :: x_low,y_low,x_hi,y_hi, xupper,yupper
    integer, allocatable :: minlevel(:,:), maxlevel(:,:)

    allocate(minlevel(mx,my), maxlevel(mx,my))
    
    minlevel = 0
    maxlevel = 0

    xupper = xlower + mx*dx
    yupper = ylower + my*dy

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

        if (ylower >= regions(m)%y_hi .or. yupper <= regions(m)%y_low) then
            cycle rloop  ! no intersection
        else
            j1 = max(floor((regions(m)%y_low - ylower) / dy) + 1, 1)
            j2 = min(floor((regions(m)%y_hi - ylower) / dy) + 1, my)
        endif

        do j=j1,j2
            do i=i1,i2
                minlevel(i,j) = max(minlevel(i,j), regions(m)%min_level)
                maxlevel(i,j) = max(maxlevel(i,j), regions(m)%max_level)
            enddo
         enddo
    enddo rloop

    do j=1,my
        do i=1,mx
         if (minlevel(i,j) > maxlevel(i,j)) then
              write(6,*) '*** Error: this should never happen!'
              write(6,*) '*** minlevel > maxlevel in flagregions'
              stop
         endif

         if (maxlevel(i,j) /= 0) then
             ! this point lies in at least one region, so may need
             ! to modify the exisiting flag at this point...
             if (level < minlevel(i,j)) then
                 ! Require refinement of this cell:
                 amrflags(i,j) = DOFLAG
             else if (level >= maxlevel(i,j)) then
                 ! Do not refine of this cell:
                 amrflags(i,j) = DONTFLAG
             ! else leave amrflags(i,j) alone.
             endif
         endif

        enddo 
    enddo 

end subroutine flagregions2
