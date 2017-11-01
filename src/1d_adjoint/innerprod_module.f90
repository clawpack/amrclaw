module innerprod_module

contains

    function calculate_innerproduct(t,q,r,mx_f,xlower_f,dx_f) result(innerprod)

        use adjoint_module

        implicit none

        real(kind=8), intent(in) :: t, xlower_f, dx_f
        integer, intent(in) :: r, mx_f
        real(kind=8), allocatable :: q_innerprod1(:), q_innerprod2(:), q_innerprod(:)
        double precision, allocatable :: q_interp(:,:), innerprod(:)
        logical, allocatable :: mask_qforward(:), mask_qadjoint(:)
        double precision, intent(in) :: q(:,:)
        integer :: mx_a, mitot_a, mptr_a
        real(kind=8) :: dx_a, xlower_a, xupper_a, xupper_f
        integer :: i, i1, i2, level, loc, x1, x2, z

        allocate(q_interp(adjoints(1)%meqn,mx_f))
        allocate(mask_qforward(mx_f))
        allocate(innerprod(mx_f))
        allocate(q_innerprod1(mx_f))
        allocate(q_innerprod2(mx_f))
        allocate(q_innerprod(mx_f))

        xupper_f = xlower_f + mx_f*dx_f
        innerprod = 0.0

        !write(*,*) "In calculate_innerproduct"
        !write(*,*) "Forward grid boundaries: ", xlower_f, xupper_f

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

            !write(*,*) " "
            !write(*,*) "Adjoint patch number ",z
            !write(*,*) "Adjoint patch boundaries: ", xlower_a, xupper_a, mx_a

            loc = adjoints(r)%loc(mptr_a)

            ! Total number of points in x
            mitot_a = mx_a + 2*adjoints(r)%nghost

            ! Check if adjoint patch overlaps with forward patch
            x1 = max(xlower_f,xlower_a)
            x2 = min(xupper_f,xupper_a)
            !write(*,*) "x1, x2: ", x1, x2

            if (x1 > x2) then
                ! Skipping interpolation if grids don't overlap
                mask_qforward = .false.
                continue
            else
                allocate(mask_qadjoint(mx_a))

                ! Create a mask that is .true. only in part of patch intersecting forward patch:
                i1 = max(int((x1 - xlower_a + 0.5d0*dx_a) / dx_a), 0)
                i2 = min(int((x2 - xlower_a + 0.5d0*dx_a) / dx_a) + 1, mx_a+1)

                if (.true.) then
                x1 = xlower_a + (i1-0.5d0)*dx_a
                x2 = xlower_a + (i2-0.5d0)*dx_a
                !write(*,*) 'patch intersecting fpatch: i1,i2: ',i1,i2,x1,x2
                endif

                forall (i=1:mx_a)
                    mask_qadjoint(i) = ((i >= i1) .and. (i <= i2))
                end forall

                ! Create a mask that is .true. only in part of forward patch intersecting patch:

                i1 = max(int((x1 - xlower_f + 0.5d0*dx_f) / dx_f), 0)
                i2 = min(int((x2 - xlower_f + 0.5d0*dx_f) / dx_f) + 1, mx_f+1)

                if (.true.) then
                x1 = xlower_f + (i1-0.5d0)*dx_f
                x2 = xlower_f + (i2-0.5d0)*dx_f
                !write(*,*) 'fpatch intersecting patch: i1,i2: ',i1,i2,x1,x2
                endif

                do i=1,mx_f
                    mask_qforward(i) = ((i >= i1) .and. (i <= i2))
                enddo

                ! Interpolate adjoint values to q_interp
                ! Note that values in q_interp will only be set properly where 
                ! mask_qadjoint == .true.
                !write(*,*) "Inner_mod, about to call interp"
                call interp_adjoint( &
                        adjoints(r)%meqn, r, q_interp, &
                        xlower_a, dx_a, mx_a, xlower_f, dx_f, mx_f, &
                        mask_qadjoint, mptr_a)

                q_innerprod1 = 0.d0
                ! For each over lapping point, calculate inner product
                forall(i = 1:mx_f, mask_qforward(i))
                    q_innerprod1(:) = abs(dot_product(q(:,i),q_interp(:,i)))
                end forall

                q_innerprod2 = 0.d0
                q_innerprod = q_innerprod1

                if (r .ne. 1) then
                    call interp_adjoint( &
                        adjoints(r)%meqn, r-1, q_interp, &
                        xlower_a, dx_a, mx_a, xlower_f, dx_f, mx_f, &
                        mask_qadjoint, mptr_a)

                    ! For each over lapping point, calculate inner product
                    forall(i = 1:mx_f, mask_qforward(i))
                        q_innerprod2(:) = abs(dot_product(q(:,i),q_interp(:,i)))
                    end forall

                    ! Assign max value to q_innerprod
                    do i=1,mx_f
                        q_innerprod(i) = max(q_innerprod1(i), q_innerprod2(i))
                    enddo
                endif

                do i=1,mx_f
                    if (q_innerprod(i) > innerprod(i)) then
                        innerprod(i) = q_innerprod(i)
                    endif
                enddo

                deallocate(mask_qadjoint)
            endif


        enddo

    end function calculate_innerproduct

end module innerprod_module
