module innerprod_module

contains

    subroutine calculate_innerproduct(q,k,mx_f,my_f,xlower_f, &
               ylower_f,dx_f,dy_f,meqn_f,mbc_f, innerprod)

        use adjoint_module

        implicit none

        real(kind=8), intent(in) :: xlower_f,ylower_f,dx_f,dy_f
        integer :: k,mx_f,my_f,meqn_f,mbc_f
        real(kind=8), intent(in) :: q(meqn_f,1-mbc_f:mx_f+mbc_f,1-mbc_f:my_f+mbc_f)

        integer :: mx_a, my_a, mptr_a, mbc_a
        integer :: i, j, i1, i2, j1, j2, level, loc, z
        real(kind=8) :: dx_a, xlower_a, xupper_a, xupper_f
        real(kind=8) :: dy_a, ylower_a, yupper_a, yupper_f
        real(kind=8) :: x1, x2, y1, y2

        real(kind=8), intent(inout) :: innerprod(mx_f,my_f)
        real(kind=8) :: q_innerprod(mx_f,my_f)
        logical :: mask_forward(mx_f,my_f)
        real(kind=8) :: q_interp(meqn_f,mx_f,my_f)

        logical, allocatable :: mask_adjoint(:,:)

        xupper_f = xlower_f + mx_f*dx_f
        yupper_f = ylower_f + my_f*dy_f

        ! Loop over patches in adjoint solution
        do z = 1, adjoints(k)%ngrids
            mptr_a = adjoints(k)%gridpointer(z)
            level = adjoints(k)%gridlevel(mptr_a)

            ! Number of points in x and y
            mx_a = adjoints(k)%ncellsx(mptr_a)
            my_a = adjoints(k)%ncellsy(mptr_a)

            ! Finding extreem values for grid
            xlower_a = adjoints(k)%xlowvals(mptr_a)
            dx_a = adjoints(k)%hxposs(level)
            xupper_a = xlower_a + mx_a*dx_a
            ylower_a = adjoints(k)%ylowvals(mptr_a)
            dy_a = adjoints(k)%hyposs(level)
            yupper_a = ylower_a + my_a*dy_a

            loc = adjoints(k)%loc(mptr_a)

            ! Check if adjoint patch overlaps with forward patch
            x1 = max(xlower_f,xlower_a)
            x2 = min(xupper_f,xupper_a)
            y1 = max(ylower_f,ylower_a)
            y2 = min(yupper_f,yupper_a)

            if ((x1 > x2) .or. (y1 > y2)) then
                ! Skipping interpolation if grids don't overlap
                mask_forward = .false.
                continue
            else
                mbc_a = adjoints(k)%nghost
                allocate(mask_adjoint(1-mbc_a:mx_a+mbc_a, 1-mbc_a:my_a+mbc_a))

                ! Create a mask that is .true. only in part of patch intersecting forward patch:
                i1 = max(int((x1 - xlower_a + 0.5d0*dx_a) / dx_a), 0)
                i2 = min(int((x2 - xlower_a + 0.5d0*dx_a) / dx_a) + 1, mx_a+1)
                j1 = max(int((y1 - ylower_a + 0.5d0*dy_a) / dy_a), 0)
                j2 = min(int((y2 - ylower_a + 0.5d0*dy_a) / dy_a) + 1, my_a+1)

                forall (i=1-mbc_a:mx_a+mbc_a, j=1-mbc_a:my_a+mbc_a)
                    mask_adjoint(i,j) = ((i >= i1) .and. (i <= i2) .and. &
                                       (j >= j1) .and. (j <= j2))
                end forall

                ! Create a mask that is .true. only in part of forward patch intersecting patch:

                i1 = max(int((x1 - xlower_f + 0.5d0*dx_f) / dx_f)+1, 0)
                i2 = min(int((x2 - xlower_f + 0.5d0*dx_f) / dx_f), mx_f)
                j1 = max(int((y1 - ylower_f + 0.5d0*dy_f) / dy_f)+1, 0)
                j2 = min(int((y2 - ylower_f + 0.5d0*dy_f) / dy_f), my_f)

                forall (i=1:mx_f, j=1:my_f)
                    mask_forward(i,j) = ((i >= i1) .and. (i <= i2) .and. &
                                       (j >= j1) .and. (j <= j2))
                end forall

                ! Interpolate adjoint values to q_interp
                ! Note that values in q_interp will only be set properly where
                ! mask_adjoint == .true.
                call interp_adjoint( &
                        adjoints(k)%meqn, k, q_interp, xlower_a, ylower_a, &
                        dx_a, dy_a, mx_a, my_a, xlower_f, ylower_f, &
                        dx_f, dy_f, mx_f, my_f, &
                        mask_adjoint, mptr_a, mask_forward)

                q_innerprod = 0.d0
                ! For each overlapping point, calculate inner product
                forall(i = 1:mx_f, j = 1:my_f, mask_forward(i,j))
                    q_innerprod(i,j) = abs(dot_product(q(:,i,j),q_interp(:,i,j)))
                end forall

                do i=1,mx_f
                    do j=1,my_f
                        if (q_innerprod(i,j) > innerprod(i,j)) then
                            innerprod(i,j) = q_innerprod(i,j)
                        endif
                    enddo
                enddo

                deallocate(mask_adjoint)
            endif
        enddo

    end subroutine calculate_innerproduct

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

end module innerprod_module
