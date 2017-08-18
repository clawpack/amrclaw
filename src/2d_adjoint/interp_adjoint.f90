! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine interp_adjoint(nvar, x, y, q, k)

        use adjoint_module, only: adjoints

        implicit none

        ! Function arguments
        integer, intent(in) :: k
        integer :: ii_c, jj_c, ii_a, jj_a, nvar, z
        real(kind=8), intent(in) :: x,y
        integer :: nx,ny,loc,locaux,level,mptr,mitot,mjtot, &
                ivar,iaux,i,j, iadd, iaddaux, iaddqeta
        real(kind=8) :: xlow, ylow, xhi, yhi, dx, dy,xm,ym, &
                x_side, x_main, y_main, y_side
        real(kind=8) :: q(nvar+1), aux_a, aux_c,q_temp1(nvar+1), &
                q_temp2(nvar+1), denom, aux_interp
        logical :: y_interp, yc_interp, ya_interp

        iadd(ivar,i,j)  = loc + ivar - 1 + adjoints(k)%meqn*((j-1)*mitot+i-1)

        do ivar=1,nvar
            q(ivar) = 0.0
        enddo

        do z = 1, adjoints(k)%ngrids
            mptr = adjoints(k)%gridpointer(z)
            level = adjoints(k)%gridlevel(mptr)

            ! Number of points in x and y (nx by ny grid)
            nx = adjoints(k)%ncellsx(mptr)
            ny = adjoints(k)%ncellsy(mptr)

            ! Finding x and y extreem values for grid
            xlow = adjoints(k)%xlowvals(mptr)
            ylow = adjoints(k)%ylowvals(mptr)
            dx = adjoints(k)%hxposs(level)
            dy = adjoints(k)%hyposs(level)
            xhi = xlow + nx*dx
            yhi = ylow + ny*dy

            loc     = adjoints(k)%loc(mptr)

            ! Total number of points in x and y
            mitot = nx + 2*adjoints(k)%nghost
            mjtot = ny + 2*adjoints(k)%nghost

            if ((x < xlow) .or. (x > xhi) .or. (y < ylow) .or. (y > yhi)) then
                ! Skipping interpolation if the point of interest
                ! is not in the current grid
                continue
            else
                xm = xlow - (adjoints(k)%nghost+0.5d0)*dx
                ym = ylow - (adjoints(k)%nghost+0.5d0)*dy

                ! Finding current cell in x (i) and y (j)
                ii_c = int((x-xm)/dx)
                jj_c = int((y-ym)/dy)

                ! Finding correct cell to interpolate with
                jj_a = int(((y-ym)/dy) + 0.5d0)
                ii_a = int(((x-xm)/dx) + 0.5d0)

                if (jj_c == jj_a .and. jj_a /= 0) then
                    jj_a = jj_a - 1
                endif
                if (ii_c == ii_a .and. ii_a /= 0) then
                    ii_a = ii_a - 1
                endif
                if (jj_a >= ny) then
                    jj_a = jj_a - 1
                endif
                if (ii_a >= nx) then
                    ii_a = ii_a - 1
                endif

                ! Interpolating in y
                y_main = ym + (jj_c + 0.5d0)*dy
                if (jj_c /= jj_a) then
                    y_interp = .true.
                    y_side = ym + (jj_a + 0.5d0)*dy
                    denom = y_side - y_main

                    do ivar=1,nvar
                        q_temp1(ivar) = &
                            ((y_side - y)/denom)*adjoints(k)%alloc(iadd(ivar,ii_c,jj_c)) &
                            + ((y - y_main)/denom)*adjoints(k)%alloc(iadd(ivar,ii_c,jj_a))
                    enddo

                    do ivar=1,nvar
                        q_temp2(ivar) = &
                            ((y_side - y)/denom)*adjoints(k)%alloc(iadd(ivar,ii_a,jj_c)) &
                            + ((y - y_main)/denom)*adjoints(k)%alloc(iadd(ivar,ii_a,jj_a))
                    enddo
                else
                    y_interp = .false.
                endif

                ! Interpolating in x
                x_main = xm + (ii_c + 0.5d0)*dx
                if (ii_c /= ii_a) then
                    x_side = xm + (ii_a + 0.5d0)*dx
                    denom = x_side - x_main

                    if(y_interp) then
                        q = ((x_side - x)/denom)*q_temp1 + ((x - x_main)/denom)*q_temp2
                    else
                        do ivar=1,nvar
                            q(ivar) = &
                                ((x_side - x)/denom)*adjoints(k)%alloc(iadd(ivar,ii_c,jj_c)) &
                                + ((x - x_main)/denom)*adjoints(k)%alloc(iadd(ivar,ii_a,jj_c))
                        enddo
                    endif
                else
                    if (y_interp) then
                        q = q_temp1
                    else
                        do ivar=1,nvar
                            q(ivar) = adjoints(k)%alloc(iadd(ivar,ii_c,jj_c))
                        enddo
                    endif
                endif
            endif

        enddo

      end subroutine interp_adjoint
