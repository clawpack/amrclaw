! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine interp_adjoint(nvar, r, q_interp, xlower_a, dx_a, &
            mx_a, xlower_f, dx_f, mx_f, mask_adjoint, mptr_a, mask_forward)

        use adjoint_module, only: adjoints

        implicit none

        ! Function arguments
        integer, intent(in) :: r, nvar, mx_a
        logical, intent(in) :: mask_adjoint(1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost)
        real(kind=8), intent(in) :: xlower_a, xlower_f
        integer, intent(in) :: mx_f, mptr_a
        real(kind=8), intent(in) :: dx_f, dx_a

        integer :: ii_c, ii_a, z,level, ik(mx_f), &
                      ivar,i, iadd, iaddaux, loc
        real(kind=8) :: q_interp(nvar,mx_f), denom
        real(kind=8) :: x_side, x_main, xm, x, xhigh_a,x_f
        real(kind=8) :: dxk(mx_f), a
        logical :: mask_forward(mx_f)

        iadd(ivar,i)  = loc + ivar - 1 + adjoints(r)%meqn*(i-1)

        q_interp = 0.0
        xhigh_a  = xlower_a + mx_a*dx_a
        loc    = adjoints(r)%loc(mptr_a)
        xm = xlower_a - (adjoints(r)%nghost+0.5d0)*dx_a

        do z = 1, mx_f
            x = xlower_f + (z - 0.5d0)*dx_f

            !write(*,*) "Considering current f point ", z, x

            ! Finding current cell in adjoint x (i)
            ii_c = int((x-xm)/dx_a)
            !write(*,*) "Overlaps with adjoint cell ", ii_c

            if (ii_c >= 1-adjoints(r)%nghost .and. ii_c <= mx_a+adjoints(r)%nghost) then
            ! Interpolate only if this cell is overlapping with grid
            if (mask_adjoint(ii_c)) then

                ! Finding correct adjoint cell to interpolate with
                ii_a = int(((x-xm)/dx_a) + 0.5d0)

                if (ii_c == ii_a .and. ii_a /= 0) then
                    ii_a = ii_a - 1
                endif
                if (ii_a >= mx_a) then
                    ii_a = ii_a - 1
                endif

                ! Interpolating in x
                x_main = xm + (ii_c - 0.5d0)*dx_a
                if (ii_c /= ii_a) then
                    x_side = xm + (ii_a - 0.5d0)*dx_a
                    denom = x_side - x_main

                    do ivar=1,nvar
                        q_interp(ivar,z) = &
                            ((x_side - x)/denom)*adjoints(r)%alloc(iadd(ivar,ii_c)) &
                            + ((x - x_main)/denom)*adjoints(r)%alloc(iadd(ivar,ii_a))
                    enddo
                else
                    do ivar=1,nvar
                        q_interp(ivar,z) = adjoints(r)%alloc(iadd(ivar,ii_c))
                    enddo
                endif

                write(*,*) "Value of q_interp: ", q_interp(:,z)
            endif
            endif
        enddo

      end subroutine interp_adjoint
