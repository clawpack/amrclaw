! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine interp_adjoint(lst, lend, nvar, x, q, k)

        use adjoint_module, only: adjoints
        use amr_reload_module

        ! Function arguments
        integer, intent(in) :: lst, lend, k
        integer :: ii_c, ii_a, nvar
        real(kind=8), intent(in) :: x
        integer :: nx,loc,level,mptr,mitot, &
                ivar,i, iadd, iaddaux
        real(kind=8) :: xlow, xhi, dx,xm, &
                x_side, x_main
        real(kind=8) :: q(nvar+1),q_temp1(nvar+1), &
                q_temp2(nvar+1), denom

        iadd(ivar,i)  = loc + ivar - 1 + nvar*(i-1)

        do ivar=1,nvar+1
            q(ivar) = 0.0
        enddo

        level = lst
65      if (level .gt. lend) go to 90
            mptr = adjoints(k)%lstart(level)
70          if (mptr .eq. 0) go to 80
                ! Number of points in x (nx)
                nx = adjoints(k)%node(ndihi,mptr) - adjoints(k)%node(ndilo,mptr) + 1

                ! Finding x extreem values for grid
                xlow = adjoints(k)%rnode(cornxlo,mptr)
                dx = adjoints(k)%hxposs(level)
                xhi = xlow + nx*dx

                loc     = adjoints(k)%node(store1, mptr)

                ! Total number of points in x
                mitot = nx + 2*nghost

                if ((x < xlow) .or. (x > xhi)) then
                    ! Skipping interpolation if the point of interest
                    ! is not in the current grid
                    continue
                else
                    xm = xlow - (nghost+0.5d0)*dx

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

            mptr = adjoints(k)%node(levelptr, mptr)
            go to 70
80      level = level + 1
        go to 65
90   continue

      end subroutine interp_adjoint
