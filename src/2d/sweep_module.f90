module sweep_module
    implicit none
    contains

subroutine x_sweep_1st_order(q, fm, fp, s_x, wave_x, meqn, mwaves, mbc, mx, my, dtdx, cflgrid) 

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, mwaves
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: s_x(mwaves, 1-mbc:mx + mbc, 2-mbc:my+mbc-1)
    real(kind=8), intent(inout) :: wave_x(meqn, mwaves, 1-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(in) :: dtdx

    ! Local variables for the Riemann solver
    integer :: i,j
    real(kind=8) :: delta1, delta2, a1, a2
    integer :: m, mw, mu, mv
    real(kind=8) :: rho, bulk, cc, zz
    real(kind=8) :: amdq(meqn), apdq(meqn)

    common /cparam/ rho,bulk,cc,zz


    ! ============================================================================
    ! Perform X-Sweeps
    do j = 0,my+1
        do i = 2-mbc, mx+mbc
            ! solve Riemann problem between cell (i-1,j) and (i,j)
            mu = 2
            mv = 3
            delta1 = q(1,i,j) - q(1,i-1,j)
            delta2 = q(mu,i,j) - q(mu,i-1,j)
            a1 = (-delta1 + zz*delta2) / (2.d0*zz)
            a2 = (delta1 + zz*delta2) / (2.d0*zz)
            !        # Compute the waves.
            wave_x(1,1,i,j) = -a1*zz
            wave_x(mu,1,i,j) = a1
            wave_x(mv,1,i,j) = 0.d0
            s_x(1,i,j) = -cc

            wave_x(1,2,i,j) = a2*zz
            wave_x(mu,2,i,j) = a2
            wave_x(mv,2,i,j) = 0.d0
            s_x(2,i,j) = cc
            do m = 1,meqn
                amdq(m) = s_x(1,i,j)*wave_x(m,1,i,j)
                apdq(m) = s_x(2,i,j)*wave_x(m,2,i,j)
                if (i >= 1 .and. i<=(mx+1)) then
                    fm(m,i,j) = fm(m,i,j) + amdq(m)
                    fp(m,i,j) = fp(m,i,j) - apdq(m)
                endif
            enddo
            ! if (mcapa > 0)  then
            !     dtdxl = dtdx / aux(mcapa,i-1,j)
            !     dtdxr = dtdx / aux(mcapa,i,j)
            ! else
            !     dtdxl = dtdx
            !     dtdxr = dtdx
            ! endif
            do mw=1,mwaves
                if (i >= 1 .and. i<=(mx+1)) then
                    ! cflgrid = dmax1(cflgrid, dtdxr*s_x(mw,i,j),-dtdxl*s_x(mw,i,j))
                    cflgrid = dmax1(cflgrid, dtdx*s_x(mw,i,j),-dtdx*s_x(mw,i,j))
                endif
            enddo
        enddo
    enddo
end subroutine x_sweep_1st_order

        subroutine x_sweep_2nd_order
        end subroutine x_sweep_2nd_order

        subroutine y_sweep_1st_order
        end subroutine y_sweep_1st_order

        subroutine y_sweep_2nd_order
        end subroutine y_sweep_2nd_order
end module sweep_module
