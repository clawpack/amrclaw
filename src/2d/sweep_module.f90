module sweep_module
    implicit none
    contains

subroutine x_sweep_1st_order(q, fm, fp, s_x, wave_x, meqn, mwaves, mbc, mx, my, dtdx, cflgrid) 

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, mwaves
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
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

subroutine x_sweep_2nd_order(qold, fm, fp, gm, gp, s_x, wave_x, meqn, mwaves, mbc, mx, my, dtdx)
end subroutine x_sweep_2nd_order

subroutine y_sweep_1st_order(q, gm, gp, s_y, wave_y, meqn, mwaves, mbc, mx, my, dtdy, cflgrid)

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, mwaves
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8) :: s_y(mwaves, 2-mbc:mx+mbc-1, 1-mbc:my + mbc)
    real(kind=8) :: wave_y(meqn, mwaves, 2-mbc:mx+mbc-1, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(in) :: dtdy

    ! Local variables for the Riemann solver
    integer :: i,j
    real(kind=8) :: delta1, delta2, a1, a2
    integer :: m, mw, mu, mv
    real(kind=8) :: rho, bulk, cc, zz
    real(kind=8) :: bmdq(meqn), bpdq(meqn)

    common /cparam/ rho,bulk,cc,zz
    ! ============================================================================
    !  y-sweeps    
    do i = 0,mx+1
        do j = 2-mbc, my+mbc
            ! solve Riemann problem between cell (i,j-1) and (i,j)
            mu = 3
            mv = 2
            delta1 = q(1,i,j) - q(1,i,j-1)
            delta2 = q(mu,i,j) - q(mu,i,j-1)
            a1 = (-delta1 + zz*delta2) / (2.d0*zz)
            a2 = (delta1 + zz*delta2) / (2.d0*zz)
            !        # Compute the waves.
            wave_y(1,1,i,j) = -a1*zz
            wave_y(mu,1,i,j) = a1
            wave_y(mv,1,i,j) = 0.d0
            s_y(1,i,j) = -cc

            wave_y(1,2,i,j) = a2*zz
            wave_y(mu,2,i,j) = a2
            wave_y(mv,2,i,j) = 0.d0
            s_y(2,i,j) = cc
            do m = 1,meqn
                bmdq(m) = s_y(1,i,j)*wave_y(m,1,i,j)
                bpdq(m) = s_y(2,i,j)*wave_y(m,2,i,j)
                gm(m,i,j) = gm(m,i,j) + bmdq(m)
                gp(m,i,j) = gp(m,i,j) - bpdq(m)
            enddo
            ! if (mcapa > 0)  then
            !     dtdyl = dtdy / aux(mcapa,i-1,j)
            !     dtdyr = dtdy / aux(mcapa,i,j)
            ! else
            !     dtdyl = dtdy
            !     dtdyr = dtdy
            ! endif
            do mw=1,mwaves
                ! cflgrid = dmax1(cflgrid, dtdyr*s_y(mw,i,j),-dtdyl*s_y(mw,i,j))
                cflgrid = dmax1(cflgrid, dtdy*s_y(mw,i,j),-dtdy*s_y(mw,i,j))
            enddo
        enddo
    enddo
end subroutine y_sweep_1st_order

subroutine y_sweep_2nd_order
end subroutine y_sweep_2nd_order
end module sweep_module
