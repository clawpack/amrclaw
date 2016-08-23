! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
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
subroutine flag2refine(mx,my,mz,mbc,meqn,maux,xlower,ylower,zlower,    &
                       dx,dy,dz,t,level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)

    use amr_module, only: lfine

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,mz,mbc,meqn,maux,level
    real(kind=8), intent(in) :: xlower,ylower,zlower,dx,dy,dz,t,tolsp
    
    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    
    ! Flagging
    real(kind=8),intent(inout) :: amrflags(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    real(kind=8), intent(in) :: DONTFLAG
    real(kind=8), intent(in) :: DOFLAG
    
    logical :: allowflag
    external allowflag

    ! Locals
    integer :: i,j,k,m
    real(kind=8) :: x_c,y_c,z_c,x_low,y_low,z_low,x_hi,y_hi,z_hi
    real(kind=8) :: dqi(meqn), dqj(meqn), dqk(meqn),dq(meqn)

    ! Initialize flags
    amrflags = DONTFLAG
    
    ! Loop over interior points on this grid
    ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
    ! This information is not needed for the default flagging based on
    ! undivided differences, but might be needed in a user's version.
    ! Note that if you want to refine only in certain space-time regions,
    ! it may be easiest to use the "regions" feature.  The flags set here or
    ! in the Richardson error estimator are potentially modified by the
    ! min_level and max_level specified in any regions.

    z_loop: do k=1,mz
        z_low = zlower + (k - 1) * dz
        z_c = zlower + (k - 0.5d0) * dz
        z_hi = zlower + k * dz

      y_loop: do j=1,my
        y_low = ylower + (j - 1) * dy
        y_c = ylower + (j - 0.5d0) * dy
        y_hi = ylower + j * dy
        
        x_loop: do i = 1,mx
            x_low = xlower + (i - 1) * dx
            x_c = xlower + (i - 0.5d0) * dx
            x_hi = xlower + i * dx

            ! -----------------------------------------------------------------
                dq = 0.d0
                dqi = abs(q(:,i+1,j,k) - q(:,i-1,j,k))
                dqj = abs(q(:,i,j+1,k) - q(:,i,j-1,k))
                dqk = abs(q(:,i,j,k+1) - q(:,i,j,k-1))
                dq = max(dq,dqi,dqj,dqk)

                ! default checks all components of undivided difference:
                !do m=1,meqn
                !    if (dq(m) > tolsp) then
                     if (abs(x_c-25.) .lt. dx) then
                        amrflags(i,j,k) = DOFLAG
                        cycle x_loop
                     endif
                !enddo

        enddo x_loop
      enddo y_loop
    enddo z_loop
    !amrflags(mx/2,my/2,mz/2) = DOFLAG
    !if (level .eq. lfine .and. t .gt. 1.) amrflags(mx/2,my/2,mz/2) = DOFLAG
   
end subroutine flag2refine
