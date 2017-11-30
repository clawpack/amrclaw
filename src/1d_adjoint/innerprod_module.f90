module innerprod_module

contains

    function calculate_innerproduct(t,q,r,mx_f,xlower_f,dx_f,meqn_f,mbc_f,lcheck) result(innerprod)

        use adjoint_module

        implicit none

        real(kind=8), intent(in) :: t, xlower_f, dx_f
        integer, intent(in) :: r, mx_f, meqn_f, mbc_f, lcheck
        real(kind=8), intent(in) :: q(meqn_f,1-mbc_f:mx_f+mbc_f)

        integer :: mx_a, mptr_a
        integer :: i, i1, i2, level, loc, z
        real(kind=8) :: dx_a, xlower_a, xupper_a, xupper_f, x1, x2

        real(kind=8) :: innerprod(1-mbc_f:mx_f+mbc_f)
        real(kind=8) :: q_innerprod(1-mbc_f:mx_f+mbc_f)
        logical :: mask_forward(1-mbc_f:mx_f+mbc_f)
        real(kind=8) :: q_interp(adjoints(r)%meqn,1-mbc_f:mx_f+mbc_f)

        logical, allocatable :: mask_adjoint(:)


        xupper_f = xlower_f + mx_f*dx_f
        innerprod = 0.0

        ! Loop over patches in adjoint solution
        do z = 1, adjoints(r)%ngrids
            mptr_a = adjoints(r)%gridpointer(z)
            level = adjoints(r)%gridlevel(mptr_a)

            ! Number of points in x (nx)
            mx_a = adjoints(r)%ncellsx(mptr_a)

            ! Finding x extreem values for grid
            xlower_a = adjoints(r)%xlowvals(mptr_a)
            dx_a = adjoints(r)%hxposs(level)
            xupper_a = xlower_a + mx_a*dx_a

            loc = adjoints(r)%loc(mptr_a)

            ! Check if adjoint patch overlaps with forward patch
            x1 = max(xlower_f,xlower_a)
            x2 = min(xupper_f,xupper_a)

            if (x1 > x2) then
                ! Skipping interpolation if grids don't overlap
                mask_forward = .false.
                continue
            else

                allocate(mask_adjoint(1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost))

                ! Create a mask that is .true. only in part of patch intersecting forward patch:
                i1 = max(int((x1 - xlower_a + 0.5d0*dx_a) / dx_a), 0)
                i2 = min(int((x2 - xlower_a + 0.5d0*dx_a) / dx_a) + 1, mx_a+1)

                forall (i=1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost)
                    mask_adjoint(i) = ((i >= i1) .and. (i <= i2))
                end forall

                ! Create a mask that is .true. only in part of forward patch intersecting patch:

                i1 = max(int((x1 - xlower_f + 0.5d0*dx_f) / dx_f)+1, 0)
                i2 = min(int((x2 - xlower_f + 0.5d0*dx_f) / dx_f), mx_f)

                do i=1-mbc_f,mx_f+mbc_f
                    mask_forward(i) = ((i >= i1) .and. (i <= i2))
                enddo

                ! Interpolate adjoint values to q_interp
                ! Note that values in q_interp will only be set properly where 
                ! mask_adjoint == .true.
                call interp_adjoint( &
                        adjoints(r)%meqn, r, q_interp, &
                        xlower_a, dx_a, mx_a, xlower_f, dx_f, mx_f, &
                        mask_adjoint, mptr_a, mask_forward, mbc_f)

                q_innerprod = 0.d0
                ! For each overlapping point, calculate inner product
                forall(i = 1-mbc_f:mx_f+mbc_f, mask_forward(i))
                    q_innerprod(i) = abs(dot_product(q(:,i),q_interp(:,i)))
                end forall

                do i=1-mbc_f,mx_f+mbc_f
                    if (q_innerprod(i) > innerprod(i)) then
                        innerprod(i) = q_innerprod(i)
                    endif
                enddo

                deallocate(mask_adjoint)
            endif
        enddo

    end function calculate_innerproduct

end module innerprod_module
