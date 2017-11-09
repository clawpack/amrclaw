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

        integer :: z,level, iz(mx_f), &
                      ivar,i, iadd, iaddaux, loc
        real(kind=8) :: q_interp(nvar,mx_f), denom
        real(kind=8) :: x, xhigh_a,x_f
        real(kind=8) :: dxz(mx_f), a
        logical :: mask_forward(mx_f)

        iadd(ivar,i)  = loc + ivar - 1 + adjoints(r)%meqn*(i-1)

        q_interp = 0.0
        xhigh_a  = xlower_a + mx_a*dx_a
        loc    = adjoints(r)%loc(mptr_a)

        do z=1,mx_f
            x = xlower_f + (z - 0.5d0)*dx_f

            iz(z) = int((x - xlower_a + 0.5d0*dx_a) / dx_a)
            dxz(z) = x - (xlower_a + (iz(z)-0.5d0)*dx_a)
        enddo

        do z = 1, mx_f
            x = xlower_f + (z - 0.5d0)*dx_f

            if (mask_forward(z)) then
            ! Interpolate only if this cell is overlapping with grid
                if (mask_adjoint(iz(z))) then

                    do ivar=1,nvar
                        a = (adjoints(r)%alloc(iadd(ivar,iz(z)+1)) &
                            - adjoints(r)%alloc(iadd(ivar,iz(z)))) / dx_a
                        q_interp(ivar,z) = &
                            adjoints(r)%alloc(iadd(ivar,iz(z))) + a*dxz(z)
                    enddo
                endif
            endif
        enddo

      end subroutine interp_adjoint
