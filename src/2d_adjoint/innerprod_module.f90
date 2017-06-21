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

        allocate(q_interp(nvar+1))

        max_innerprod = 0.d0
        ! Select adjoint data
        aloop: do r=1,totnum_adjoints

          if (r .ne. 1) then
              t_nm = adjoints(r-1)%time
          else
              t_nm = 0.d0
          endif

          if ((t+adjoints(r)%time) >= trange_start .and. &
              (t+adjoints(r)%time) <=trange_final) then

            q_innerprod1 = 0.d0
            q_innerprod2 = 0.d0

            call interp_adjoint(1, adjoints(r)%lfine, nvar, &
                 x_c,y_c,q_interp,r)
            q_innerprod1 = abs( q1 * q_interp(1) &
                  + q2 * q_interp(2) + q3 * q_interp(3))

            if (r .ne. 1) then
                call interp_adjoint(1, adjoints(r)%lfine, nvar, &
                     x_c,y_c,q_interp, r-1)
                q_innerprod2 = abs(q1 * q_interp(1) &
                    + q2 * q_interp(2) + q3 * q_interp(3))
            endif

            ! Assign max value to q_innerprod
            q_innerprod = max(q_innerprod1, q_innerprod2)
            if (q_innerprod > max_innerprod) then
                max_innerprod = q_innerprod
            endif

          endif
        enddo aloop

    end function calculate_max_innerproduct

end module innerprod_module
