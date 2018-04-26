! ::::::::::::::::::::: flagregions ::::::::::::::::::::::::::::::::::
!
!> Modify array of flagged points to respect minlevels and maxlevels
!! specified by regions.
!
!
! amrflags  = array to be flagged with either the value
!             DONTFLAG (no refinement needed)  or
!             DOFLAG   (refinement desired)    
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine flagregions(mx,my,mbuff,xlower,ylower,dx,dy,level,t, &
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
    
    logical :: allowflag
    external allowflag

    ! Locals
    integer :: i,j,m,minlevel,maxlevel
    real(kind=8) :: x_low,y_low,x_hi,y_hi

    
    ! Loop over interior points on this grid
    ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi]
    do j=1,my
        y_low = ylower + (j - 1) * dy
        y_hi = ylower + j * dy
        
        do i = 1,mx
            x_low = xlower + (i - 1) * dx
            x_hi = xlower + i * dx

            minlevel = 0
            maxlevel = 0

            do m=1,num_regions
                if (t >= regions(m)%t_low .and. t <= regions(m)%t_hi .and. &
                   x_hi > regions(m)%x_low .and. x_low < regions(m)%x_hi .and. &
                   y_hi > regions(m)%y_low .and. y_low < regions(m)%y_hi ) then
                        minlevel = max(minlevel, regions(m)%min_level)
                        maxlevel = max(maxlevel, regions(m)%max_level)
                endif
             enddo

             if (minlevel > maxlevel) then
                  write(6,*) '*** Error: this should never happen!'
                  write(6,*) '*** minlevel > maxlevel in flagregions'
                  stop
             endif

             if (maxlevel /= 0) then
                 ! this point lies in at least one region, so may need
                 ! to modify the exisiting flag at this point...
                 if (level < minlevel) then
                     ! Require refinement of this cell:
                     amrflags(i,j) = DOFLAG
                 else if (level >= maxlevel) then
                     ! Do not refine of this cell:
                     amrflags(i,j) = DONTFLAG
                 ! else leave amrflags(i,j) alone.
                 endif
              endif

        enddo 
    enddo 
   
end subroutine flagregions
