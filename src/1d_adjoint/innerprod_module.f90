module innerprod_module

contains

    function calculate_max_innerproduct(t,x_c,q1,q2) result(max_innerprod)

        use adjoint_module

        real(kind=8), intent(in) :: t
        integer :: r
        real(kind=8) :: q_innerprod1, q_innerprod2, q_innerprod, max_innerprod
        double precision, allocatable :: q_interp(:)
        real(kind=8) :: x_c,q1,q2
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
            q_interp = interpolate_adjoint(1,adjoints(r)%lfine,nvar,x_c,r)
            q_innerprod1 = abs( q1 * q_interp(1) &
                  + q2 * q_interp(2))

            q_innerprod2 = 0.d0
            q_innerprod = q_innerprod1
            if (r .ne. 1) then
                q_interp = interpolate_adjoint(1,adjoints(r)%lfine,nvar,x_c,r-1)

                q_innerprod2 = abs(q1 * q_interp(1) &
                    + q2 * q_interp(2))

                ! Assign max value to q_innerprod
                q_innerprod = max(q_innerprod1, q_innerprod2)
            endif

            if (q_innerprod > max_innerprod) then
                max_innerprod = q_innerprod
            endif

          endif
        enddo aloop

    end function calculate_max_innerproduct

    function interpolate_adjoint(lst, lend, nvar, x, k) result(q)

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

    end function interpolate_adjoint

end module innerprod_module
