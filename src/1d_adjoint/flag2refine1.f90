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

subroutine flag2refine1(mx,mbc,mbuff,meqn,maux,xlower,dx,t,level, &
                            tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)

    use innerprod_module, only: calculate_max_innerproduct
    use adjoint_module, only: innerprod_index
    use regions_module

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,dx,t,tolsp
    
    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc)
    
    ! Flagging
    real(kind=8),intent(inout) :: amrflags(1-mbuff:mx+mbuff)
    real(kind=8), intent(in) :: DONTFLAG
    real(kind=8), intent(in) :: DOFLAG
    
    logical :: allowflag
    external allowflag

    ! Locals
    integer :: i,m
    real(kind=8) :: x_c,x_low,x_hi
    real(kind=8) :: dqi(meqn), dq(meqn)

    ! Initialize flags
    amrflags = DONTFLAG

      x_loop: do i = 1,mx
        x_low = xlower + (i - 1) * dx
        x_c = xlower + (i - 0.5d0) * dx
        x_hi = xlower + i * dx

        aux(innerprod_index,i) =  &
           calculate_max_innerproduct(t,x_c,q(1,i),q(2,i))

        if (aux(innerprod_index,i) > tolsp) then
            amrflags(i) = DOFLAG
            cycle x_loop
        endif

      enddo x_loop

end subroutine flag2refine1
