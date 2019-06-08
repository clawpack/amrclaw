! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
!> \callgraph
!! \callergraph
!! User routine to control flagging of points for refinement.
!!
!! Default version computes spatial difference dq in each direction and
!! for each component of q and flags any point where this is greater than
!! the tolerance tolsp.  
!! This is consistent with what the routine errsp did in
!!  earlier versions of amrclaw (4.2 and before).
!!
!! This routine can be copied to an application directory and modified to
!! implement some other desired refinement criterion.
!!
!! Points may also be flagged for refining based on a Richardson estimate
!! of the error, obtained by comparing solutions on the current grid and a
!! coarsened grid.  Points are flagged if the estimated error is larger than
!! the parameter tol in amr2ez.data, provided flag_richardson is .true.,
!! otherwise the coarsening and Richardson estimation is not performed!  
!! Points are flagged via Richardson in a separate routine.
!!
!! Once points are flagged via this routine and/or Richardson, the subroutine
!! flagregions is applied to check each point against the min_level and
!! max_level of refinement specified in any "region" set by the user.
!! So flags set here might be over-ruled by region constraints.
!!
!! **output**: amrflags
!!
!! \param mx number of cells in *i* direction
!! \param my number of cells in *j* direction
!! \param mbc width of ghost cell region
!! \param mbuff width of buffer region
!! \param meqn number of equations for the system
!! \param maux number of auxiliary variables
!! \param xlower x-coordinate of left physical boundary
!! \param ylower y-coordinate of lower physical boundary
!! \param dx spacing in *i* direction
!! \param dy spacing in *j* direction
!! \param t simulation time on this grid
!! \param level AMR level of this grid
!! \param tolsp tolerance specified by user in input file amr.data, used in default
!!         version of this routine as a tolerance for spatial differences
!! \param q grid values including ghost cells (bndry vals at specified
!!          time have already been set, so can use ghost cell values too)
!! \param aux auxiliary array on this grid patch
!! \param amrflags array to be flagged with either the value **DONTFLAG** or **DOFLAG** for each cell. 
!!        It is enlarged from grid size to include buffer regions around the grid. 
!! \param DONTFLAG value to be assigned to amrflags for cells that need no refinement
!! \param DOFLAG value to be assigned to amrflags for cells that do need refinement
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine flag2refine2(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,dx,dy,t,level, &
                            tolsp,q,aux,amrflags)

    use regions_module
    use adjoint_module, only: totnum_adjoints,innerprod_index, &
                              adjoint_flagging,select_snapshots
    use adjointsup_module, only: calculate_innerproduct
    use amr_module, only : DOFLAG, UNSET

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp
    
    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Flagging
    real(kind=8),intent(inout) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)
    
    logical :: allowflag
    external allowflag

    ! Locals
    integer :: i,j,m,r
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    real(kind=8) :: dqi(meqn), dqj(meqn), dq(meqn)

    ! Adjoint method specific variables
    logical mask_selecta(totnum_adjoints)

    ! Don't initialize flags, since they were already 
    ! flagged by flagregions2
    ! amrflags = DONTFLAG

    if(adjoint_flagging) then
        aux(innerprod_index,:,:) = 0.0
        call select_snapshots(t,mask_selecta)

        ! Loop over adjoint snapshots
        aloop: do r=1,totnum_adjoints

            ! Consider only snapshots that are within the desired time range
            if (mask_selecta(r)) then
                ! Calculate inner product with current snapshot
                call calculate_innerproduct(q,r,mx,my,xlower,   &
                        ylower,dx,dy,meqn,mbc,maux,aux)
            endif

        enddo aloop
    endif
    
    ! Loop over interior points on this grid
    ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
    ! This information is not needed for the default flagging based on
    ! undivided differences, but might be needed in a user's version.
    ! Note that if you want to refine only in certain space-time regions,
    ! it may be easiest to use the "regions" feature. The flags set here or
    ! in the Richardson error estimator are modifing the flags set by
    ! min_level and max_level specified in any regions.

    y_loop: do j=1,my
        !y_low = ylower + (j - 1) * dy
        !y_c = ylower + (j - 0.5d0) * dy
        !y_hi = ylower + j * dy
        
        x_loop: do i = 1,mx
            !x_low = xlower + (i - 1) * dx
            !x_c = xlower + (i - 0.5d0) * dx
            !x_hi = xlower + i * dx

            ! -----------------------------------------------------------------
            ! Only check undivided differences if flag hasn't been set yet. 
            ! If flag == DONTFLAG then refinement is forbidden by a region, 
            ! if flag == DOFLAG checking is not needed
            if(amrflags(i,j) == UNSET) then

                if(adjoint_flagging) then
                    if (aux(innerprod_index,i,j) > tolsp) then
                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                else
                    dq = 0.d0
                    dqi = abs(q(:,i+1,j) - q(:,i-1,j))
                    dqj = abs(q(:,i,j+1) - q(:,i,j-1))
                    dq = max(dq,dqi,dqj)

                    ! default checks all components of undivided difference:
                    do m=1,meqn
                        if (dq(m) > tolsp) then
                            amrflags(i,j) = DOFLAG
                            cycle x_loop
                        endif
                    enddo
                endif

            endif

        enddo x_loop
    enddo y_loop
   
end subroutine flag2refine2
