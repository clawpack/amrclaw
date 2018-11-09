! :::::::::: bc1amr ::::::::::::::::::::::::::::::::::::::::::::::;
!
!     Take a grid patch with mesh width hx of dimensions nrow
!     and set the values of any piece of
!     of the patch which extends outside the physical domain 
!     using the boundary conditions. 
!
!     ------------------------------------------------
!     # Standard boundary condition choices for amr in clawpack
!
!     # At each boundary  k = 1 (left),  2 (right):
!     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
!     #            =  1  for zero-order extrapolation
!     #            =  2  for periodic boundary conditions
!     #            =  3  for solid walls, assuming this can be implemented
!     #                  by reflecting the data about the boundary and then
!     #                  negating the 2'nd (for k=1,2) component of q.
!     #            =  5  sphere bcs (left half maps to right half of same 
!     #                  side, and vice versa), as if domain folded in half
!     ------------------------------------------------
!
!     The edges of the grid patch are at
!        xlo_patch  -- left edge
!        xhi_patch --  right edge
!
!     The physical domain itself is a rectangle bounded by
!        xlower  -- left edge
!        xupper  -- right edge
!     
!     the picture is the following: 
!            (xlower)                    (xupper)
!        |______|_________|_________________|
!        |      |         |                 |
!   (xlo_patch)       (xhi_patch)
!        
!
!     Any cells that lie outside the physical domain are ghost cells whose
!     values should be set in this routine.  This is tested for by comparing
!     xlo_patch with xlower to see if values need to be set at the left, as in
!     the figure above, and similarly at the other boundary.
!
!     Patches are guaranteed to have at least 1 row of cells filled
!     with interior values so it is possible to  extrapolate. 
!     Fix trimbd if you want more than 1 row pre-set.
!
!     Make sure the order the boundaries are specified is correct
!     so that diagonal corner cells are also properly taken care of.
!
!     Periodic boundaries are set before calling this routine, so if the
!     domain is periodic in one direction only you
!     can safely extrapolate in the other direction. 
!
!     Don't overwrite ghost cells in periodic directions!
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

subroutine bc1amr(val, aux, nrow, meqn,naux, hx, level, time,           &
                  xlo_patch, xhi_patch) 

    use amr_module, only: mthbc, xlower, xupper
    use amr_module, only: xperdom

    implicit none

    ! Input/Output
    integer, intent(in) :: nrow, meqn, naux, level
    real(kind=8), intent(in) :: hx, time
    real(kind=8), intent(in) :: xlo_patch, xhi_patch
    real(kind=8), intent(in out) :: val(meqn, nrow)
    real(kind=8), intent(in out) :: aux(naux, nrow)
    
    ! Local storage
    integer :: i, ibeg, nxl, nxr
    real(kind=8) :: hxmarg

    hxmarg = hx * .01d0

    ! Use periodic boundary condition specialized code only, if only one 
    ! boundary is periodic we still proceed below
    if (xperdom) then
        return
    end if

    ! Each check has an initial check to ensure that the boundary is a real
    ! boundary condition and otherwise skips the code.  Otherwise 
    !-------------------------------------------------------
    ! Left boundary:
    !-------------------------------------------------------
    if (xlo_patch < xlower-hxmarg) then
        ! number of grid cells from this patch lying outside physical domain:
        nxl = int((xlower + hxmarg - xlo_patch) / hx)

        select case(mthbc(1))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            case(1) ! Zero-order extrapolation
                do i=1, nxl
                    val(:, i) = val(:, nxl + 1)
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do i=1, nxl
                    aux(:, i) = aux(:, 2 * nxl + 1 - i)
                    val(:, i) = val(:, 2 * nxl + 1 - i)
                end do
                ! negate the normal velocity:
                do i=1, nxl
                    val(2, i) = -val(2, i)
                end do

            case default
                print *, "Invalid boundary condition requested."
                stop
        end select
    end if

    !-------------------------------------------------------
    ! Right boundary:
    !-------------------------------------------------------
    if (xhi_patch > xupper+hxmarg) then

        ! number of grid cells lying outside physical domain:
        nxr = int((xhi_patch - xupper + hxmarg) / hx)
        ibeg = max(nrow - nxr + 1, 1)

        select case(mthbc(2))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            case(1) ! Zero-order extrapolation
                do i = ibeg, nrow
                    val(:, i) = val(:, ibeg - 1)
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do i=ibeg, nrow
                    val(:, i) = val(:, 2 * ibeg - 1 - i)
                end do
                ! negate the normal velocity:
                do i = ibeg, nrow
                    val(2, i) = -val(2, i)
                end do

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

end subroutine bc1amr