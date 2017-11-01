! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine interp_adjoint(nvar, r, q_interp, xlower_a, dx_a, &
            mx_a, xlower_f, dx_f, mx_f, mask_qadjoint, mptr_a)

        use adjoint_module, only: adjoints

        implicit none

        ! Function arguments
        integer, intent(in) :: r, nvar
        logical, intent(in) :: mask_qadjoint(mx_a)
        real(kind=8), intent(in) :: xlower_a, xlower_f
        integer, intent(in) :: mx_f, mx_a
        real(kind=8), intent(in) :: dx_f, dx_a

        integer :: ii_c, ii_a, z,level,mptr_a,mitot_a, &
                      ivar,i, iadd, iaddaux, loc
        real(kind=8) :: q_interp(nvar,mx_f), denom
        real(kind=8) :: x_side, x_main, xm, x, xhigh_a

        iadd(ivar,i)  = loc + ivar - 1 + adjoints(r)%meqn*(i-1)

        q_interp = 0.0
        xhigh_a  = xlower_a + mx_a*dx_a
        loc    = adjoints(r)%loc(mptr_a)
        ! Total number of points in x
        mitot_a  = mx_a + 2*adjoints(r)%nghost
        xm = xlower_a !- (adjoints(r)%nghost+0.5d0)*dx_a

        !write(*,*) "In interp_adjoint"

        do z = 1, mx_f
            x = xlower_f + (z - 0.5d0)*dx_f

            !write(*,*) "Considering current f point ", z, x

            ! Finding current cell in adjoint x (i)
            ii_c = int((x-xm)/dx_a)
            !write(*,*) "Overlaps with adjoint cell ", ii_c

            ! Interpolate only if this cell is overlapping with grid
            if (mask_qadjoint(ii_c)) then

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

        enddo

      end subroutine interp_adjoint
