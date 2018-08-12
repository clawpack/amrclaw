! :::::::::: BC3AMR ::::::::::::::::::::::::::::::::::::::::::::::;
!
!     Take a grid patch with mesh widths hx,hy,hz, of dimensions nrow by
!     ncol by nfil,  and set the values of any piece of
!     of the patch which extends outside the physical domain
!     using the boundary conditions.
!     ------------------------------------------------
!     # Standard boundary condition choices for amr3ez in clawpack
!
!     # At each boundary  k = 1 (xlower),  2 (xupper),  3 (ylower ), 4 (yupper),
!                         5 (zlower) 6 (zupper)
!     #  
!     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
!     #            =  1  for zero-order extrapolation
!     #            =  2  for periodic boundary coniditions
!     #            =  3  for solid walls, assuming this can be implemented
!     #                  by reflecting the data about the boundary and then
!     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
!     #                  or 4'th (for k = 5,6) component of q.
!     ------------------------------------------------
!
!     The corners of the grid patch are at
!        (xlo_patch,ylo_patch,zlo_patch)  --  lower left corner
!        (xhi_patch,yhi_patch,zhi_patch)  --  upper right corner
!
!     The physical domain itself is a rectangular parallelopiped bounded by
!        (xlower,ylower,zlower)  -- lower front left corner
!        (xupper,yupper,zupper)  -- upper rear right corner
!
!     the picture is the following:
!
!                            __________________________(xupper,yupper,zupper)
!                           /                         /|
!                          /                         / |
!                         /                         /  |
!                        /_________________________/   |
!                        |                         |   |
!                        |                         |   |
!                     ___|_____(xhi_patch,yhi_patch,zhi_patch) 
!                    /___|____/|                   |   |
!                    |   |    ||                   |   |
!                    |   |    ||                   |   |
!                    |   |    ||                   |   |
!                    |___|____|/                   |   |
!  (xlo_patch,ylo_patch,zlo_patch)                 |  /                       
!                        |                         | /
!                        |_________________________|/
!  (xlower,ylower,zlower)
!
!     Any cells that lie outside the physical domain are ghost cells whose
!     values should be set in this routine.  This is tested for by comparing
!     xlo_patch with xlower to see if values need to be set at the left, as in
!     the figure above, and similarly at the other boundaries.
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

subroutine bc3amr(val,aux,nrow,ncol,nfil,meqn,naux, hx, hy, hz, level, time,  &
                  xlo_patch, xhi_patch,                     &
                  ylo_patch, yhi_patch,                     &
                  zlo_patch, zhi_patch) 

    use amr_module, only: mthbc, xlower, ylower, xupper, yupper, zlower, zupper
    use amr_module, only: xperdom,yperdom,zperdom

    implicit none

    ! Input/Output
    integer, intent(in) :: nrow, ncol, nfil, meqn, naux, level
    real(kind=8), intent(in) :: hx, hy, hz, time
    real(kind=8), intent(in) :: xlo_patch, xhi_patch
    real(kind=8), intent(in) :: ylo_patch, yhi_patch
    real(kind=8), intent(in) :: zlo_patch, zhi_patch
    real(kind=8), intent(in out) :: val(meqn, nrow, ncol, nfil)
    real(kind=8), intent(in out) :: aux(naux, nrow, ncol, nfil)
    
    ! Local storage
    integer :: i, j, k, ibeg, jbeg, kbeg, nxl, nxr, nyl, nyr, nzl, nzr
    real(kind=8) :: hxmarg, hymarg, hzmarg

    hxmarg = hx * .01d0
    hymarg = hy * .01d0
    hzmarg = hz * .01d0

    ! Use periodic boundary condition specialized code only, if only one 
    ! boundary is periodic we still proceed below
    if (xperdom .and. yperdom .and. zperdom) then
        return
    end if

    ! Each check has an initial check to ensure that the boundary is a real
    ! boundary condition and otherwise skips the code.  Otherwise 
    !-------------------------------------------------------
    ! xlower boundary:
    !-------------------------------------------------------
    if (xlo_patch < xlower-hxmarg) then
        ! number of grid cells from this patch lying outside physical domain:
        nxl = int((xlower + hxmarg - xlo_patch) / hx)

        select case(mthbc(1))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            case(1) ! Zero-order extrapolation
                do k = 1, nfil
                    do j = 1, ncol
                        do i=1, nxl
                            val(:, i, j, k) = val(:, nxl + 1, j, k)
                        end do
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do k = 1, nfil
                    do j = 1, ncol
                        do i=1, nxl
                            val(:, i, j, k) = val(:, 2 * nxl + 1 - i, j, k)
                        end do
                    end do
                end do
                ! negate the normal velocity:
                do k = 1, nfil
                    do j = 1, ncol
                        do i=1, nxl
                            val(2, i, j, k) = -val(2, i, j, k)
                        end do
                    end do
                end do

            case default
                print *, "Invalid boundary condition requested."
                stop
        end select
    end if

    !-------------------------------------------------------
    ! xupper boundary:
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
                do k = 1, nfil
                    do i = ibeg, nrow
                        do j = 1, ncol
                            val(:, i, j, k) = val(:, ibeg - 1, j, k)
                        end do
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do k = 1, nfil
                    do i=ibeg, nrow
                        do j = 1, ncol
                            val(:, i, j, k) = val(:, 2 * ibeg - 1 - i, j, k)
                        end do
                    end do
                end do
                ! negate the normal velocity:
                do k = 1, nfil
                    do i = ibeg, nrow
                        do j = 1, ncol
                            val(2, i, j, k) = -val(2, i, j, k)
                        end do
                    end do
                end do

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! ylower boundary:
    !-------------------------------------------------------
    if (ylo_patch < ylower - hymarg) then

        ! number of grid cells lying outside physical domain:
        nyl = int((ylower + hymarg - ylo_patch) / hy)

        select case(mthbc(3))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."
            
            case(1) ! Zero-order extrapolation
                do k = 1, nfil
                    do j = 1, nyl
                        do i = 1, nrow
                            val(:, i ,j, k) = val(:, i, nyl + 1, k)
                        end do
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do k = 1, nfil
                    do j = 1, nyl
                        do i = 1, nrow
                            val(:, i ,j, k) = val(:, i, 2 * nyl + 1 - j, k)
                        end do
                    end do
                end do
                ! negate the normal velocity:
                do k = 1, nfil
                    do j = 1, nyl
                        do i = 1, nrow
                            val(3, i ,j, k) = -val(3, i, j, k)
                        end do
                    end do
                end do

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! yupper boundary:
    !-------------------------------------------------------
    if (yhi_patch > yupper + hymarg) then

        ! number of grid cells lying outside physical domain:
        nyr = int((yhi_patch - yupper + hymarg) / hy)
        jbeg = max(ncol - nyr + 1, 1)

        select case(mthbc(4))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."

            case(1) ! Zero-order extrapolation
                do k = 1, nfil
                    do j = jbeg, ncol
                        do i = 1, nrow
                            val(:, i, j, k) = val(:, i, jbeg - 1, k)
                        end do
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do k = 1, nfil
                    do j = jbeg, ncol 
                        do i = 1, nrow
                            val(:, i, j, k) = val(:, i, 2 * jbeg - 1 - j, k)
                        end do
                    end do
                end do
                ! negate the normal velocity:
                do k = 1, nfil
                    do j = jbeg, ncol
                        do i = 1, nrow
                            val(3, i, j, k) = -val(3, i, j, k)
                        end do
                    end do
                end do

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! Z lower boundary:
    !-------------------------------------------------------
    if (zlo_patch < zlower - hzmarg) then

        ! Numnber of ghost cells lying outside physical domain
        nzl = (zlower + hzmarg - zlo_patch) / hz

        select case(mthbc(5))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."

            case(1) ! Zero-order extrapolation
                do k = 1, nzl
                    do j = 1, ncol
                        do i = 1, nrow
                            val(:, i, j, k) = val(:, i, j, nzl + 1)
                        end do
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do k = 1, nzl
                    do j = 1, ncol 
                        do i = 1, nrow
                            val(:, i, j, k) = val(:, i, j, 2 * nzl + 1 - k)
                        end do
                    end do
                end do
                ! negate the normal velocity:
                do k = 1, nzl
                    do j = 1, ncol
                        do i = 1, nrow
                            val(4, i, j, k) = -val(4, i, j, k)
                        end do
                    end do
                end do

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

    !-------------------------------------------------------
    ! Z upper boundary:
    !-------------------------------------------------------
    if (zhi_patch > zupper + hzmarg) then

        ! Numnber of ghost cells lying outside physical domain
        nzr = (zhi_patch - zupper + hzmarg) / hz
        kbeg = max(nfil - nzr + 1, 1)

        select case(mthbc(6))
            case(0) ! User defined boundary condition
                ! Replace this code with a user defined boundary condition
                stop "A user defined boundary condition was not provided."

            case(1) ! Zero-order extrapolation
                do k = kbeg, nfil
                    do j = 1, ncol
                        do i = 1, nrow
                            val(:, i, j, k) = val(:, i, j, kbeg - 1)
                        end do
                    end do
                end do

            case(2) ! Periodic boundary condition
                continue

            case(3) ! Wall boundary conditions
                do k = kbeg, nfil
                    do j = 1, ncol 
                        do i = 1, nrow
                            val(:, i, j, k) = val(:, i, j, 2 * kbeg - 1 - k)
                        end do
                    end do
                end do
                ! negate the normal velocity:
                do k = kbeg, nfil
                    do j = 1, ncol
                        do i = 1, nrow
                            val(4, i, j, k) = -val(4, i, j, k)
                        end do
                    end do
                end do

            case default
                print *, "Invalid boundary condition requested."
                stop

        end select
    end if

end subroutine bc3amr