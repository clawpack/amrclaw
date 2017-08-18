! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine interp_adjoint(nvar, x, q, k)

        use adjoint_module, only: adjoints

        implicit none

        ! Function arguments
        integer, intent(in) :: k
        integer :: ii_c, ii_a, nvar, z
        real(kind=8), intent(in) :: x
        integer :: nx,loc,level,mptr,mitot, &
                ivar,i, iadd, iaddaux
        real(kind=8) :: xlow, xhi, dx,xm, &
                x_side, x_main
        real(kind=8) :: q(nvar+1),q_temp1(nvar+1), &
                q_temp2(nvar+1), denom

        iadd(ivar,i)  = loc + ivar - 1 + adjoints(k)%meqn*(i-1)

        do ivar=1,nvar
            q(ivar) = 0.0
        enddo

        do z = 1, adjoints(k)%ngrids
            mptr = adjoints(k)%gridpointer(z)
            level = adjoints(k)%gridlevel(mptr)

            ! Number of points in x (nx)
            nx = adjoints(k)%ncellsx(mptr)

            ! Finding x extreem values for grid
            xlow = adjoints(k)%xlowvals(mptr)
            dx = adjoints(k)%hxposs(level)
            xhi = xlow + nx*dx

            loc     = adjoints(k)%loc(mptr)

            ! Total number of points in x
            mitot = nx + 2*adjoints(k)%nghost

            if ((x < xlow) .or. (x > xhi)) then
                ! Skipping interpolation if the point of interest
                ! is not in the current grid
                continue
            else
                xm = xlow - (adjoints(k)%nghost+0.5d0)*dx

                ! Finding current cell in x (i)
                ii_c = int((x-xm)/dx)

                ! Finding correct cell to interpolate with
                ii_a = int(((x-xm)/dx) + 0.5d0)

                if (ii_c == ii_a .and. ii_a /= 0) then
                    ii_a = ii_a - 1
                endif
                if (ii_a >= nx) then
                    ii_a = ii_a - 1
                endif


                ! Interpolating in x
                x_main = xm + (ii_c + 0.5d0)*dx
                if (ii_c /= ii_a) then
                    x_side = xm + (ii_a + 0.5d0)*dx
                    denom = x_side - x_main

                    do ivar=1,nvar
                        q(ivar) = &
                            ((x_side - x)/denom)*adjoints(k)%alloc(iadd(ivar,ii_c)) &
                            + ((x - x_main)/denom)*adjoints(k)%alloc(iadd(ivar,ii_a))
                    enddo
                else
                    do ivar=1,nvar
                        q(ivar) = adjoints(k)%alloc(iadd(ivar,ii_c))
                    enddo
                endif
            endif

        enddo

      end subroutine interp_adjoint
