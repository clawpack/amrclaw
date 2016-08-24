module innerprod_module

contains

    function calculate_max_innerproduct(t,x_c,y_c,q1,q2,q3) result(max_innerprod)

        use adjoint_module

        real(kind=8), intent(in) :: t
        integer :: r
        real(kind=8) :: q_innerprod1, q_innerprod2, q_innerprod, max_innerprod
        double precision, allocatable :: q_interp(:)
        real(kind=8) :: x_c,y_c,q1,q2,q3
        real(kind=8) :: t_nm

        max_innerprod = 0.d0
        ! Select adjoint data
        aloop: do r=1,totnum_adjoints

          if (r .ne. 1) then
              t_nm = adjoints(r-1)%time
          else
              t_nm = 0.d0
          endif

          if (t <= adjoints(r)%time .and. &
              min((t + (tfinal - trange_start)),tfinal) >= t_nm) then

            !write (*,*) "About to interpolate."
            q_interp = interpolate_adjoint(1,adjoints(r)%lfine,nvar,x_c,y_c,r)
            q_innerprod1 = abs( q1 * q_interp(1) &
                  + q2 * q_interp(2) + q3 * q_interp(3))

            q_innerprod2 = 0.d0
            q_innerprod = q_innerprod1
            if (r .ne. 1) then
                q_interp = interpolate_adjoint(1,adjoints(r)%lfine,nvar,x_c,y_c,r-1)

                q_innerprod2 = abs(q1 * q_interp(1) &
                    + q2 * q_interp(2) + q3 * q_interp(3))

                ! Assign max value to q_innerprod
                q_innerprod = max(q_innerprod1, q_innerprod2)
            endif

            if (q_innerprod > max_innerprod) then
                max_innerprod = q_innerprod
            endif

          endif
        enddo aloop

    end function calculate_max_innerproduct

    function interpolate_adjoint(lst, lend, nvar, x, y, k) result(q)

        use adjoint_module, only: adjoints
        use amr_reload_module

        ! Function arguments
        integer, intent(in) :: lst, lend, k
        integer :: ii_c, jj_c, ii_a, jj_a, nvar
        real(kind=8), intent(in) :: x,y
        integer :: nx,ny,loc,locaux,level,mptr,mitot,mjtot, &
                ivar,iaux,i,j, iadd, iaddaux, iaddqeta
        real(kind=8) :: xlow, ylow, xhi, yhi, dx, dy,xm,ym, &
                x_side, x_main, y_main, y_side
        real(kind=8) :: q(nvar+1), aux_a, aux_c,q_temp1(nvar+1), &
                q_temp2(nvar+1), denom, aux_interp
        logical :: y_interp, yc_interp, ya_interp

        iadd(ivar,i,j)  = loc + ivar - 1 + nvar*((j-1)*mitot+i-1)

        do ivar=1,nvar+1
            q(ivar) = 0.0
        enddo

        level = lst
65      if (level .gt. lend) go to 90
            mptr = adjoints(k)%lstart(level)
70          if (mptr .eq. 0) go to 80
                ! Number of points in x and y (nx by ny grid)
                nx = adjoints(k)%node(ndihi,mptr) - adjoints(k)%node(ndilo,mptr) + 1
                ny = adjoints(k)%node(ndjhi,mptr) - adjoints(k)%node(ndjlo,mptr) + 1

                ! Finding x and y extreem values for grid
                xlow = adjoints(k)%rnode(cornxlo,mptr)
                ylow = adjoints(k)%rnode(cornylo,mptr)
                dx = adjoints(k)%hxposs(level)
                dy = adjoints(k)%hyposs(level)
                xhi = xlow + nx*dx
                yhi = ylow + ny*dy

                loc     = adjoints(k)%node(store1, mptr)

                ! Total number of points in x and y
                mitot = nx + 2*nghost
                mjtot = ny + 2*nghost

                if ((x < xlow) .or. (x > xhi) .or. (y < ylow) .or. (y > yhi)) then
                    ! Skipping interpolation if the point of interest
                    ! is not in the current grid
                    continue
                else
                    xm = xlow - (nghost+0.5d0)*dx
                    ym = ylow - (nghost+0.5d0)*dy

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

            mptr = adjoints(k)%node(levelptr, mptr)
            go to 70
80      level = level + 1
        go to 65
90   continue

    end function interpolate_adjoint

end module innerprod_module
