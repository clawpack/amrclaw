! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine interp_adjoint(nvar, r, q_interp, xlower_a, ylower_a, dx_a, dy_a, &
            mx_a, my_a, xlower_f, ylower_f, dx_f, dy_f, mx_f, my_f, &
            mask_adjoint, mptr_a, mask_forward)

        use adjoint_module, only: adjoints

        implicit none

        ! Function arguments
        integer, intent(in) :: r, nvar, mx_a, my_a
        logical, intent(in) :: mask_adjoint(1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost, &
                                       1-adjoints(r)%nghost:my_a+adjoints(r)%nghost)
        real(kind=8), intent(in) :: xlower_a, xlower_f, ylower_a, ylower_f
        integer, intent(in) :: mx_f, my_f, mptr_a
        real(kind=8), intent(in) :: dx_f, dx_a, dy_f, dy_a

        integer :: z,k, iz, jk, mitot
        integer :: ivar, i, j, iadd, iaddaux, loc
        real(kind=8) :: q_interp(nvar,mx_f,my_f), denom
        real(kind=8) :: x, xhigh_a, y, yhigh_a
        real(kind=8) :: dxz,dyk, a, b, c
        logical :: mask_forward(mx_f,my_f)

        iadd(ivar,i,j)  = loc + ivar - 1 + adjoints(r)%meqn*((j-1)*mitot+i-1)

        q_interp = 0.0
        xhigh_a  = xlower_a + mx_a*dx_a
        yhigh_a = ylower_a + my_a*dx_a
        loc    = adjoints(r)%loc(mptr_a)
        mitot = adjoints(r)%ncellsx(mptr_a) + 2*adjoints(r)%nghost

        do z = 1,mx_f
          do k = 1,my_f
            if (mask_forward(z,k)) then
                x = xlower_f + (z - 0.5d0)*dx_f
                y = ylower_f + (k - 0.5d0)*dy_f

                !TODO: Why does iz and jk have an added 1 at the end?

                iz = int((x - xlower_a + 0.5d0*dx_a) / dx_a) + 1
                dxz = x - (xlower_a + (iz-0.5d0)*dx_a)
                jk = int((y - ylower_a + 0.5d0*dy_a) / dy_a) + 1
                dyk = y - (ylower_a + (jk-0.5d0)*dy_a)

                ! Interpolate only if this cell is overlapping with grid
                if (mask_adjoint(iz,jk)) then
                    do ivar=1,nvar

                        a = (adjoints(r)%alloc(iadd(ivar,iz+1,jk)) &
                             - adjoints(r)%alloc(iadd(ivar,iz,jk))) / dx_a
                        b = (adjoints(r)%alloc(iadd(ivar,iz,jk+1)) &
                             - adjoints(r)%alloc(iadd(ivar,iz,jk))) / dy_a
                        c = (adjoints(r)%alloc(iadd(ivar,iz+1,jk+1)) &
                             + adjoints(r)%alloc(iadd(ivar,iz,jk)) &
                             - adjoints(r)%alloc(iadd(ivar,iz+1,jk)) &
                             - adjoints(r)%alloc(iadd(ivar,iz,jk+1))) / (dx_a*dy_a)

                        q_interp(ivar,z,k) = &
                             adjoints(r)%alloc(iadd(ivar,iz,jk)) &
                             + a*dxz + b*dyk + c*dxz*dyk

                    enddo
                endif
            endif
          enddo
        enddo

      end subroutine interp_adjoint
