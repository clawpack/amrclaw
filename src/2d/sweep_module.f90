#include "amr_macros.H"

module sweep_module

    use amr_module
    use problem_para_module, only: rho,bulk,cc,zz
    use cuda_module, only: max_reduce_device_2d, max_reduce_device_local_2d

    implicit none

    contains


subroutine x_sweep_1st_order(q, fm, fp, s_x, wave_x, meqn, mwaves, mbc, mx, my, dtdx, cflgrid) 

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, mwaves
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: s_x(mwaves, 2-mbc:mx + mbc, 2-mbc:my+mbc-1)
    real(kind=8), intent(inout) :: wave_x(meqn, mwaves, 2-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(in) :: dtdx

    ! Local variables for the Riemann solver
    integer :: i,j
    real(kind=8) :: delta1, delta2, a1, a2
    integer :: m, mw, mu, mv
    real(kind=8) :: amdq(meqn), apdq(meqn)


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

attributes(global) &
subroutine x_sweep_1st_order_gpu(q, fm, fp, s_x, wave_x, mbc, mx, my, dtdx, cfls, ngrids, id, cc, zz) 

    implicit none

    integer, value, intent(in) :: mbc, mx, my
    real(kind=8), intent(in) ::     q(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: fm(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: fp(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) ::    s_x(2-mbc:mx + mbc, 2-mbc:my+mbc-1, NWAVES)
    real(kind=8), intent(inout) :: wave_x(2-mbc:mx+mbc, 2-mbc:my+mbc-1, NEQNS, NWAVES)
    real(kind=8), value, intent(in) :: dtdx
    real(kind=8), value, intent(in) :: cc, zz
    integer, value, intent(in) :: ngrids, id
    real(kind=8), intent(inout) :: cfls(ngrids,2)

    ! Local variables for the Riemann solver
    integer :: i,j, tidx, tidy
    real(kind=8) :: delta1, delta2, a1, a2
    integer :: m, mw
    real(kind=8) :: amdq(NEQNS), apdq(NEQNS)
    real(kind=8) :: cfl_local
    real(kind=8) :: atomic_result


    double precision, shared :: cfl_s(blockDim%x, blockDim%y)


    tidx = threadIdx%x
    tidy = threadIdx%y
    i = (blockIdx%x-1) * blockDim%x + threadIdx%x
    j = (blockIdx%y-1) * blockDim%y + threadIdx%y
    ! we shift i and j such that they are mapped to the loop:
    ! do j = 0,my+1
    !     do i = 2-mbc, mx+mbc
    i = i + (2-mbc) - 1 ! now i = 1 is mapped to i = 2-mbc
    j = j + 0 - 1 ! now j = 1 is mapped to j = 0

    if (i > (mx+mbc) .or. j > (my+1) ) then
        return
    endif

    cfl_s(tidx, tidy) = 0.d0



    ! ============================================================================
    delta1 = q(i,j,1) - q(i-1,j, 1)
    delta2 = q(i,j,2) - q(i-1,j, 2)
    a1 = (-delta1 + zz*delta2) / (2.d0*zz)
    a2 = (delta1 + zz*delta2) / (2.d0*zz)

    !        # Compute the waves.
    s_x(i,j,1) = -cc
    s_x(i,j,2) = cc

    wave_x(i,j,1,1) = -a1*zz
    wave_x(i,j,2,1) = a1
    wave_x(i,j,3,1) = 0.d0
    wave_x(i,j,1,2) = a2*zz
    wave_x(i,j,2,2) = a2
    wave_x(i,j,3,2) = 0.d0

    amdq(:) = s_x(i,j,1)*wave_x(i,j,:,1)
    apdq(:) = s_x(i,j,2)*wave_x(i,j,:,2)
    fm(i,j,:) = fm(i,j,:) + amdq(:)
    fp(i,j,:) = fp(i,j,:) - apdq(:)

    do mw=1,NWAVES
        if (i >= 1 .and. i<=(mx+1)) then
            cfl_s(tidx, tidy) = dmax1(cfl_s(tidx,tidy), dtdx*s_x(i,j,mw),-dtdx*s_x(i,j,mw))
        endif
    enddo

    call syncthreads()
    call max_reduce_device_local_2d(cfl_s, mx+2*mbc-1, my+2, cfl_local)

    ! Write to a global cfl
    if ( tidx == 1 .and. tidy==1) then
        atomic_result = atomicmax(cfls(id,1), cfl_local)
    endif

end subroutine x_sweep_1st_order_gpu

subroutine x_sweep_2nd_order(fm, fp, gm, gp, s_x, wave_x, meqn, mwaves, mbc, mx, my, dtdx)

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, mwaves
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: s_x(mwaves, 2-mbc:mx + mbc, 2-mbc:my+mbc-1)
    real(kind=8), intent(in) :: wave_x(meqn, mwaves, 2-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8), intent(in) :: dtdx

    ! Local variables for the Riemann solver
    real(kind=8) :: wave_x_tilde(meqn, mwaves, 2-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: cqxx(meqn)
    real(kind=8) :: amdq(meqn), apdq(meqn)
    real(kind=8) :: bpamdq(meqn), bmamdq(meqn), bpapdq(meqn), bmapdq(meqn)
    real(kind=8) :: delta1, delta2, a1, a2
    real(kind=8) :: dot, wnorm2, wlimitr, abs_sign, c, r
    integer :: i,j
    integer :: m, mw, mu, mv
    logical limit


    limit = .false.
    do mw=1,mwaves
        if (mthlim(mw) .gt. 0) limit = .true.
    enddo

    if (method(2).ne.1) then ! if second-order
        wave_x_tilde = wave_x
    endif
    !     -----------------------------------------------------------
    !     # modify F fluxes for second order q_{xx} correction terms
    !     # and solve for transverse waves
    !     -----------------------------------------------------------
    do j = 0,my+1 ! my loop
        do i = 1, mx+1 ! mx loop
            if (method(2).ne.1) then ! if second-order
                ! # apply limiter to waves:
                if (limit) then ! limiter if
                    do mw=1,mwaves ! mwaves loop
                        if (mthlim(mw) .eq. 0) cycle
                        dot = 0.d0
                        wnorm2 = 0.d0
                        do m=1,meqn
                            wnorm2 = wnorm2 + wave_x(m,mw,i,j)**2
                        enddo
                        if (wnorm2.eq.0.d0) cycle

                        if (s_x(mw,i,j) .gt. 0.d0) then
                            do m=1,meqn
                                dot = dot + wave_x(m,mw,i,j)*wave_x(m,mw,i-1,j)
                            enddo
                        else
                            do m=1,meqn
                                dot = dot + wave_x(m,mw,i,j)*wave_x(m,mw,i+1,j)
                            enddo
                        endif

                        r = dot / wnorm2

                        ! choose limiter
                        if (mthlim(mw) .eq. 1) then
                            !               --------
                            !               # minmod
                            !               --------
                            wlimitr = dmax1(0.d0, dmin1(1.d0, r))

                        else if (mthlim(mw) .eq. 2) then
                            !               ----------
                            !               # superbee
                            !               ----------
                            wlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))

                        else if (mthlim(mw) .eq. 3) then
                            !               ----------
                            !               # van Leer
                            !               ----------
                            wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))

                        else if (mthlim(mw) .eq. 4) then
                            !               ------------------------------
                            !               # monotinized centered
                            !               ------------------------------
                            c = (1.d0 + r)/2.d0
                            wlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
                        else if (mthlim(mw) .eq. 5) then
                            !               ------------------------------
                            !               # Beam-Warming
                            !               ------------------------------
                            wlimitr = r
                        else
                            print *, 'Unrecognized limiter.'
                            stop
                        endif
                        !
                        !  # apply limiter to waves:
                        !
                        do m=1,meqn
                            wave_x_tilde(m,mw,i,j) = wlimitr * wave_x_tilde(m,mw,i,j)
                        enddo
                    enddo ! end mwave loop
                endif ! end limiter if

                ! if (mcapa > 0)  then
                !     dtdxl = dtdx / aux(mcapa,i-1,j)
                !     dtdxr = dtdx / aux(mcapa,i,j)
                ! else
                !     dtdxl = dtdx
                !     dtdxr = dtdx
                ! endif
                ! dtdxave = 0.5d0 * (dtdxl + dtdxr)

                ! second order corrections:
                do m=1,meqn
                    cqxx(m) = 0.d0
                    do mw=1,mwaves
                        if (use_fwaves) then
                            abs_sign = dsign(1.d0,s_x(mw,i,j))
                        else
                            abs_sign = dabs(s_x(mw,i,j))
                        endif

                        cqxx(m) = cqxx(m) + abs_sign * &
                            ! (1.d0 - dabs(s_x(mw,i,j))*dtdxave) * wave_x_tilde(m,mw,i,j)
                            (1.d0 - dabs(s_x(mw,i,j))*dtdx) * wave_x_tilde(m,mw,i,j)
                    enddo
                    fp(m,i,j) = fp(m,i,j) + 0.5d0 * cqxx(m)
                    fm(m,i,j) = fm(m,i,j) + 0.5d0 * cqxx(m)
                enddo
            endif ! end if second-order 
    ! ##### solve for transverse waves and add to gp and gm
            if (method(3).ne.0) then ! if transverse propagation
                ! reconstruct amdq and apdq
                do m=1,meqn
                    amdq(m) = s_x(1,i,j)*wave_x(m,1,i,j)
                    apdq(m) = s_x(2,i,j)*wave_x(m,2,i,j)
                enddo
                if (method(2).gt.1 .and. method(3).eq.2) then
                    ! incorporate cqxx into amdq and apdq so that it is split also.
                    do m=1,meqn
                        amdq(m) = amdq(m) + cqxx(m)
                        apdq(m) = apdq(m) - cqxx(m)
                    enddo
                endif

                mu = 2
                mv = 3
                ! ##### solve for bpamdq and bmamdq
                a1 = (-amdq(1) + zz*amdq(mv)) / (2.d0*zz)
                a2 = (amdq(1) + zz*amdq(mv)) / (2.d0*zz)
                !        # The down-going flux difference bmasdq is the product  -c * wave
                bmamdq(1) = cc * a1*zz
                bmamdq(mu) = 0.d0
                bmamdq(mv) = -cc * a1
                !        # The up-going flux difference bpasdq is the product  c * wave
                bpamdq(1) = cc * a2*zz
                bpamdq(mu) = 0.d0
                bpamdq(mv) = cc * a2

                if (mcapa > 0)  then
                    !todo
                    print * , "To be implemented"
                else
                    do m =1,meqn
                        gm(m,i-1,j) = gm(m,i-1,j) - 0.5d0*dtdx * bmamdq(m)
                        gp(m,i-1,j) = gp(m,i-1,j) - 0.5d0*dtdx * bmamdq(m)

                        gm(m,i-1,j+1) = gm(m,i-1,j+1) - 0.5d0*dtdx * bpamdq(m)
                        gp(m,i-1,j+1) = gp(m,i-1,j+1) - 0.5d0*dtdx * bpamdq(m)
                    enddo
                endif

                ! # solve for bpapdq and bmapdq
                a1 = (-apdq(1) + zz*apdq(mv)) / (2.d0*zz)
                a2 = (apdq(1) + zz*apdq(mv)) / (2.d0*zz)
                !        # The down-going flux difference bmasdq is the product  -c * wave
                bmapdq(1) = cc * a1*zz
                bmapdq(mu) = 0.d0
                bmapdq(mv) = -cc * a1
                !        # The up-going flux difference bpasdq is the product  c * wave
                bpapdq(1) = cc * a2*zz
                bpapdq(mu) = 0.d0
                bpapdq(mv) = cc * a2

                if (mcapa > 0)  then
                    !todo
                    print * , "To be implemented"
                else
                    do m =1,meqn
                        gm(m,i,j) = gm(m,i,j) - 0.5d0*dtdx * bmapdq(m)
                        gp(m,i,j) = gp(m,i,j) - 0.5d0*dtdx * bmapdq(m)

                        gm(m,i,j+1) = gm(m,i,j+1) - 0.5d0*dtdx * bpapdq(m)
                        gp(m,i,j+1) = gp(m,i,j+1) - 0.5d0*dtdx * bpapdq(m)
                    enddo
                endif
            endif ! end if transverse propagation
    ! ##### END solve for transverse waves and add to gp and gm
        enddo ! end mx loop
    enddo ! end my loop
end subroutine x_sweep_2nd_order


! For now, this kernel assumes 
! 1) limit == true
! 2) use 2nd order
! 3) include 2nd-order terms in transverse waves
! 4) van Leer limiter is used
! 5) use_fwaves == false
! 6) mcapa = 0
! We should write separate kernels for cases where 
! the options above are different

attributes(global) &
subroutine x_sweep_2nd_order_gpu(fm, fp, gm, gp, s_x, wave_x, mbc, mx, my, dtdx, cc, zz)

    implicit none

    integer, value, intent(in) :: mbc, mx, my
    real(kind=8), intent(inout) :: fm(1-mbc:mx+mbc, 1-mbc:my+mbc,NEQNS)
    real(kind=8), intent(inout) :: fp(1-mbc:mx+mbc, 1-mbc:my+mbc,NEQNS)
    real(kind=8), intent(inout) :: gm(1-mbc:mx+mbc, 1-mbc:my+mbc,NEQNS)
    real(kind=8), intent(inout) :: gp(1-mbc:mx+mbc, 1-mbc:my+mbc,NEQNS)
    real(kind=8), intent(in) ::    s_x(2-mbc:mx+mbc, 2-mbc:my+mbc-1, NWAVES)
    real(kind=8), intent(in) :: wave_x(2-mbc:mx+mbc, 2-mbc:my+mbc-1, NEQNS, NWAVES)

    real(kind=8), value, intent(in) :: dtdx
    real(kind=8), value, intent(in) :: cc, zz

    ! Local variables for the Riemann solver
    real(kind=8) :: cqxx(NEQNS)
    real(kind=8) :: amdq(NEQNS), apdq(NEQNS)
    real(kind=8) :: wave_x_tilde(NEQNS, NWAVES)
    real(kind=8) :: bpamdq(NEQNS), bmamdq(NEQNS), bpapdq(NEQNS), bmapdq(NEQNS)
    real(kind=8) :: a1, a2
    real(kind=8) :: dot, wnorm2, wlimitr, r
    integer :: i,j
    integer :: m, mw
    real(kind=8) :: atomic_result


    i = (blockIdx%x-1) * blockDim%x + threadIdx%x
    j = (blockIdx%y-1) * blockDim%y + threadIdx%y
    ! we shift i and j such that they are mapped to the loop:
    ! do j = 0,my+1 ! my loop
    !     do i = 1, mx+1 ! mx loop
    ! i = 1 already corresponds to i = 1
    j = j + 0 - 1 ! now j = 1 is mapped to j = 0

    if (i > (mx+1) .or. j > (my+1) ) then
        return
    endif
    !     -----------------------------------------------------------
    !     # modify F fluxes for second order q_{xx} correction terms
    !     # and solve for transverse waves
    !     -----------------------------------------------------------
    do mw=1,NWAVES ! mwaves loop
        dot = 0.d0
        wnorm2 = 0.d0
        do m=1,NEQNS
            wnorm2 = wnorm2 + wave_x(i,j,m,mw)**2
        enddo
        if (wnorm2.eq.0.d0) then
            wave_x_tilde(:,mw) =  wave_x(i,j,:,mw)
        else

            if (s_x(i,j,mw) .gt. 0.d0) then
                do m=1,NEQNS
                    dot = dot + wave_x(i,j,m,mw)*wave_x(i-1,j,m,mw)
                enddo
            else
                do m=1,NEQNS
                    dot = dot + wave_x(i,j,m,mw)*wave_x(i+1,j,m,mw)
                enddo
            endif

            r = dot / wnorm2
            !               ----------
            !               # van Leer
            !               ----------
            wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
            !
            !  # apply limiter to waves:
            !
            wave_x_tilde(:,mw) = wlimitr * wave_x(i,j,:,mw)
        endif
    enddo ! end mwave loop

    cqxx = 0.d0
    do mw = 1, NWAVES
        cqxx(:) = cqxx(:) + &
            dabs(s_x(i,j,mw)) * (1.d0 - dabs(s_x(i,j,mw))*dtdx) * wave_x_tilde(:,mw)
    enddo

    fp(i,j,:) = fp(i,j,:) + 0.5d0 * cqxx(:)
    fm(i,j,:) = fm(i,j,:) + 0.5d0 * cqxx(:)

    ! ##### solve for transverse waves and add to gp and gm
    ! reconstruct amdq and apdq
    amdq(:) = s_x(i,j,1)*wave_x(i,j,:,1)
    apdq(:) = s_x(i,j,2)*wave_x(i,j,:,2)
    ! incorporate cqxx into amdq and apdq so that it is split also.
    amdq(:) = amdq(:) + cqxx(:)
    apdq(:) = apdq(:) - cqxx(:)

    ! ##### solve for bpamdq and bmamdq
    a1 = (-amdq(1) + zz*amdq(3)) / (2.d0*zz)
    !        # The down-going flux difference bmasdq is the product  -c * wave
    bmamdq(1) = cc * a1*zz
    bmamdq(2) = 0.d0
    bmamdq(3) = -cc * a1

    !        # The up-going flux difference bpasdq is the product  c * wave
    a2 = (amdq(1) + zz*amdq(3)) / (2.d0*zz)
    bpamdq(1) = cc * a2*zz
    bpamdq(2) = 0.d0
    bpamdq(3) = cc * a2

    do m =1,NEQNS
        atomic_result = atomicadd(gm(i-1,j  ,m), - 0.5d0*dtdx * bmamdq(m))
        atomic_result = atomicadd(gp(i-1,j  ,m), - 0.5d0*dtdx * bmamdq(m))
        atomic_result = atomicadd(gm(i-1,j+1,m), - 0.5d0*dtdx * bpamdq(m))
        atomic_result = atomicadd(gp(i-1,j+1,m), - 0.5d0*dtdx * bpamdq(m))
    enddo

    ! # solve for bpapdq and bmapdq
    a1 = (-apdq(1) + zz*apdq(3)) / (2.d0*zz)
    !        # The down-going flux difference bmasdq is the product  -c * wave
    bmapdq(1) = cc * a1*zz
    bmapdq(2) = 0.d0
    bmapdq(3) = -cc * a1
    !        # The up-going flux difference bpasdq is the product  c * wave
    a2 = (apdq(1) + zz*apdq(3)) / (2.d0*zz)
    bpapdq(1) = cc * a2*zz
    bpapdq(2) = 0.d0
    bpapdq(3) = cc * a2

    do m =1,NEQNS
        atomic_result = atomicadd(gm(i,j  ,m), - 0.5d0*dtdx * bmapdq(m))
        atomic_result = atomicadd(gp(i,j  ,m), - 0.5d0*dtdx * bmapdq(m))
        atomic_result = atomicadd(gm(i,j+1,m), - 0.5d0*dtdx * bpapdq(m))
        atomic_result = atomicadd(gp(i,j+1,m), - 0.5d0*dtdx * bpapdq(m))
    enddo
    return
end subroutine x_sweep_2nd_order_gpu

subroutine x_sweep_2nd_order_simple(fm, fp, gm, gp, s_x, wave_x, meqn, mwaves, mbc, mx, my, dtdx)

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, mwaves
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: s_x(mwaves, 2-mbc:mx + mbc, 2-mbc:my+mbc-1)
    real(kind=8), intent(in) :: wave_x(meqn, mwaves, 2-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8), intent(in) :: dtdx

    ! Local variables for the Riemann solver
    real(kind=8) :: wave_x_tilde(meqn, mwaves, 2-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: cqxx(meqn)
    real(kind=8) :: amdq(meqn), apdq(meqn)
    real(kind=8) :: bpamdq(meqn), bmamdq(meqn), bpapdq(meqn), bmapdq(meqn)
    real(kind=8) :: delta1, delta2, a1, a2
    real(kind=8) :: dot, wnorm2, wlimitr, abs_sign, c, r
    integer :: i,j
    integer :: m, mw

    !     -----------------------------------------------------------
    !     # modify F fluxes for second order q_{xx} correction terms
    !     # and solve for transverse waves
    !     -----------------------------------------------------------
    ! wave_x_tilde = wave_x
    do j = 0,my+1 ! my loop
        do i = 1, mx+1 ! mx loop
            do mw=1,mwaves ! mwaves loop
                dot = 0.d0
                wnorm2 = 0.d0
                do m=1,meqn
                    wnorm2 = wnorm2 + wave_x(m,mw,i,j)**2
                enddo
                if (wnorm2.eq.0.d0) then
                    wave_x_tilde(:,mw,i,j) =  wave_x(:,mw,i,j)
                else

                    if (s_x(mw,i,j) .gt. 0.d0) then
                        do m=1,meqn
                            dot = dot + wave_x(m,mw,i,j)*wave_x(m,mw,i-1,j)
                        enddo
                    else
                        do m=1,meqn
                            dot = dot + wave_x(m,mw,i,j)*wave_x(m,mw,i+1,j)
                        enddo
                    endif

                    r = dot / wnorm2
                    !               ----------
                    !               # van Leer
                    !               ----------
                    wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))
                    do m=1,meqn
                        wave_x_tilde(m,mw,i,j) = wlimitr * wave_x(m,mw,i,j)
                    enddo
                endif
            enddo ! end mwave loop
                
            ! second order corrections:
            do m=1,meqn
                cqxx(m) = 0.d0
                do mw=1,mwaves
                    cqxx(m) = cqxx(m) + & 
                        dabs(s_x(mw,i,j)) * (1.d0 - dabs(s_x(mw,i,j))*dtdx) * wave_x_tilde(m,mw,i,j)
                enddo
            enddo
            do m=1,meqn
                fp(m,i,j) = fp(m,i,j) + 0.5d0 * cqxx(m)
                fm(m,i,j) = fm(m,i,j) + 0.5d0 * cqxx(m)
            enddo

    ! ##### solve for transverse waves and add to gp and gm
            ! reconstruct amdq and apdq
            do m=1,meqn
                amdq(m) = s_x(1,i,j)*wave_x(m,1,i,j)
                apdq(m) = s_x(2,i,j)*wave_x(m,2,i,j)
            enddo
            ! incorporate cqxx into amdq and apdq so that it is split also.
            do m=1,meqn
                amdq(m) = amdq(m) + cqxx(m)
                apdq(m) = apdq(m) - cqxx(m)
            enddo

            ! ##### solve for bpamdq and bmamdq
            a1 = (-amdq(1) + zz*amdq(3)) / (2.d0*zz)
            a2 = (amdq(1) + zz*amdq(3)) / (2.d0*zz)
            !        # The down-going flux difference bmasdq is the product  -c * wave
            bmamdq(1) = cc * a1*zz
            bmamdq(2) = 0.d0
            bmamdq(3) = -cc * a1
            !        # The up-going flux difference bpasdq is the product  c * wave
            bpamdq(1) = cc * a2*zz
            bpamdq(2) = 0.d0
            bpamdq(3) = cc * a2

            do m =1,meqn
                gm(m,i-1,j) = gm(m,i-1,j) - 0.5d0*dtdx * bmamdq(m)
                gp(m,i-1,j) = gp(m,i-1,j) - 0.5d0*dtdx * bmamdq(m)

                gm(m,i-1,j+1) = gm(m,i-1,j+1) - 0.5d0*dtdx * bpamdq(m)
                gp(m,i-1,j+1) = gp(m,i-1,j+1) - 0.5d0*dtdx * bpamdq(m)
            enddo

            ! # solve for bpapdq and bmapdq
            a1 = (-apdq(1) + zz*apdq(3)) / (2.d0*zz)
            a2 = (apdq(1) + zz*apdq(3)) / (2.d0*zz)
            !        # The down-going flux difference bmasdq is the product  -c * wave
            bmapdq(1) = cc * a1*zz
            bmapdq(2) = 0.d0
            bmapdq(3) = -cc * a1
            !        # The up-going flux difference bpasdq is the product  c * wave
            bpapdq(1) = cc * a2*zz
            bpapdq(2) = 0.d0
            bpapdq(3) = cc * a2

            do m =1,meqn
                gm(m,i,j) = gm(m,i,j) - 0.5d0*dtdx * bmapdq(m)
                gp(m,i,j) = gp(m,i,j) - 0.5d0*dtdx * bmapdq(m)

                gm(m,i,j+1) = gm(m,i,j+1) - 0.5d0*dtdx * bpapdq(m)
                gp(m,i,j+1) = gp(m,i,j+1) - 0.5d0*dtdx * bpapdq(m)
            enddo
        enddo ! end mx loop
    enddo ! end my loop
end subroutine x_sweep_2nd_order_simple

subroutine y_sweep_1st_order(q, gm, gp, s_y, wave_y, meqn, mwaves, mbc, mx, my, dtdy, cflgrid)

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, mwaves
    real(kind=8), intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: s_y(mwaves, 2-mbc:mx+mbc-1, 2-mbc:my + mbc)
    real(kind=8), intent(inout) :: wave_y(meqn, mwaves, 2-mbc:mx+mbc-1, 2-mbc:my+mbc)
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(in) :: dtdy

    ! Local variables for the Riemann solver
    integer :: i,j
    real(kind=8) :: delta1, delta2, a1, a2
    integer :: m, mw, mu, mv
    real(kind=8) :: bmdq(meqn), bpdq(meqn)

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


attributes(global) &
subroutine y_sweep_1st_order_gpu(q, gm, gp, s_y, wave_y, mbc, mx, my, dtdy, cfls, ngrids, id, cc, zz) 

    implicit none

    integer, value, intent(in) :: mbc, mx, my
    real(kind=8), intent(in) ::     q(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: gm(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: gp(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) ::    s_y(2-mbc:mx+mbc-1, 2-mbc:my + mbc, NWAVES)
    real(kind=8), intent(inout) :: wave_y(2-mbc:mx+mbc-1, 2-mbc:my+mbc, NEQNS, NWAVES)
    real(kind=8), value, intent(in) :: dtdy
    real(kind=8), value, intent(in) :: cc, zz
    integer, value, intent(in) :: ngrids, id
    real(kind=8), intent(inout) :: cfls(ngrids,2)

    ! Local variables for the Riemann solver
    integer :: i,j, tidx, tidy
    real(kind=8) :: delta1, delta2, a1, a2
    integer :: m, mw
    real(kind=8) :: bmdq(NEQNS), bpdq(NEQNS)
    real(kind=8) :: cfl_local
    real(kind=8) :: atomic_result

    double precision, shared :: cfl_s(blockDim%x, blockDim%y)


    tidx = threadIdx%x
    tidy = threadIdx%y
    i = (blockIdx%x-1) * blockDim%x + threadIdx%x
    j = (blockIdx%y-1) * blockDim%y + threadIdx%y
    ! we shift i and j such that they are mapped to the loop:
    ! do i = 0,mx+1
    !     do j = 2-mbc, my+mbc
    i = i + 0 - 1 ! now i = 1 is mapped to i = 2-mbc
    j = j + (2-mbc) - 1 ! now j = 1 is mapped to j = 0

    if (i > (mx+1) .or. j > (my+mbc) ) then
        return
    endif

    cfl_s(tidx, tidy) = 0.d0



    ! ============================================================================
    delta1 = q(i,j,1) - q(i,j-1,1)
    delta2 = q(i,j,3) - q(i,j-1,3)
    a1 = (-delta1 + zz*delta2) / (2.d0*zz)
    a2 = (delta1 + zz*delta2) / (2.d0*zz)
    !        # Compute the waves.
    s_y(i,j,1) = -cc
    s_y(i,j,2) = cc

    wave_y(i,j,1,1) = -a1*zz
    wave_y(i,j,2,1) = 0.d0
    wave_y(i,j,3,1) = a1

    wave_y(i,j,1,2) = a2*zz
    wave_y(i,j,2,2) = 0.d0
    wave_y(i,j,3,2) = a2

    bmdq(:) = s_y(i,j,1)*wave_y(i,j,:,1)
    bpdq(:) = s_y(i,j,2)*wave_y(i,j,:,2)
    gm(i,j,:) = gm(i,j,:) + bmdq(:)
    gp(i,j,:) = gp(i,j,:) - bpdq(:)
    do mw=1,NWAVES
        if (j >= 1 .and. j<=(my+1)) then
            cfl_s(tidx, tidy) = dmax1(cfl_s(tidx,tidy), dtdy*s_y(i,j,mw),-dtdy*s_y(i,j,mw))
        endif
    enddo

    call syncthreads()
    call max_reduce_device_local_2d(cfl_s, mx+2, my+2*mbc-1, cfl_local)

    ! Write to a global cfl
    if ( tidx == 1 .and. tidy==1) then
        atomic_result = atomicmax(cfls(id,2), cfl_local)
    endif

end subroutine y_sweep_1st_order_gpu


subroutine y_sweep_2nd_order(fm, fp, gm, gp, s_y, wave_y, meqn, mwaves, mbc, mx, my, dtdy)

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, mwaves
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(in) :: s_y(mwaves, 2-mbc:mx+mbc-1, 2-mbc:my + mbc)
    real(kind=8), intent(in) :: wave_y(meqn, mwaves, 2-mbc:mx+mbc-1, 2-mbc:my+mbc)
    real(kind=8), intent(in) :: dtdy

    ! Local variables for the Riemann solver
    real(kind=8) :: wave_y_tilde(meqn, mwaves, 2-mbc:mx+mbc-1, 2-mbc:my+mbc)
    real(kind=8) :: cqyy(meqn)
    real(kind=8) :: bmdq(meqn), bpdq(meqn)
    real(kind=8) :: apbmdq(meqn), ambmdq(meqn), apbpdq(meqn), ambpdq(meqn)
    real(kind=8) :: delta1, delta2, a1, a2
    real(kind=8) :: dot, wnorm2, wlimitr, abs_sign, c, r
    integer :: i,j
    integer :: m, mw, mu, mv
    logical limit


    limit = .false.
    do mw=1,mwaves
        if (mthlim(mw) .gt. 0) limit = .true.
    enddo

    if (method(2).ne.1) then ! if second-order
        wave_y_tilde = wave_y
    endif
    !     -----------------------------------------------------------
    !     # modify G fluxes for second order q_{yy} correction terms:
    !     # and transverse waves
    !     -----------------------------------------------------------
    do i = 0, mx+1 ! mx loop
        do j = 1, my+1 ! my loop
            if (method(2).ne.1) then ! if second-order
                ! # apply limiter to waves:
                if (limit) then ! limiter if
                    do mw=1,mwaves ! mwaves loop
                        if (mthlim(mw) .eq. 0) cycle
                        dot = 0.d0
                        wnorm2 = 0.d0
                        do m=1,meqn
                            wnorm2 = wnorm2 + wave_y(m,mw,i,j)**2
                        enddo
                        if (wnorm2.eq.0.d0) cycle

                        if (s_y(mw,i,j) .gt. 0.d0) then
                            do m=1,meqn
                                dot = dot + wave_y(m,mw,i,j)*wave_y(m,mw,i,j-1)
                            enddo
                        else
                            do m=1,meqn
                                dot = dot + wave_y(m,mw,i,j)*wave_y(m,mw,i,j+1)
                            enddo
                        endif

                        r = dot / wnorm2

                        ! choose limiter
                        if (mthlim(mw) .eq. 1) then
                            !               --------
                            !               # minmod
                            !               --------
                            wlimitr = dmax1(0.d0, dmin1(1.d0, r))

                        else if (mthlim(mw) .eq. 2) then
                            !               ----------
                            !               # superbee
                            !               ----------
                            wlimitr = dmax1(0.d0, dmin1(1.d0, 2.d0*r), dmin1(2.d0, r))

                        else if (mthlim(mw) .eq. 3) then
                            !               ----------
                            !               # van Leer
                            !               ----------
                            wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))

                        else if (mthlim(mw) .eq. 4) then
                            !               ------------------------------
                            !               # monotinized centered
                            !               ------------------------------
                            c = (1.d0 + r)/2.d0
                            wlimitr = dmax1(0.d0, dmin1(c, 2.d0, 2.d0*r))
                        else if (mthlim(mw) .eq. 5) then
                            !               ------------------------------
                            !               # Beam-Warming
                            !               ------------------------------
                            wlimitr = r
                        else
                            print *, 'Unrecognized limiter.'
                            stop
                        endif
                        !
                        !  # apply limiter to waves:
                        !
                        do m=1,meqn
                            wave_y_tilde(m,mw,i,j) = wlimitr * wave_y_tilde(m,mw,i,j)
                        enddo
                    enddo ! end mwave loop
                endif ! end limiter if
                ! if (mcapa > 0)  then
                !     dtdyl = dtdy / aux(mcapa,i,j-1)
                !     dtdyr = dtdy / aux(mcapa,i,j)
                ! else
                !     dtdyl = dtdy
                !     dtdyr = dtdy
                ! endif
                ! dtdyave = 0.5d0 * (dtdyl + dtdyr)

                ! second order corrections:
                do m=1,meqn
                    cqyy(m) = 0.d0
                    do mw=1,mwaves
                        if (use_fwaves) then
                            abs_sign = dsign(1.d0,s_y(mw,i,j))
                        else
                            abs_sign = dabs(s_y(mw,i,j))
                        endif

                        cqyy(m) = cqyy(m) + abs_sign * &
                            ! (1.d0 - dabs(s_y(mw,i,j))*dtdyave) * wave_y_tilde(m,mw,i,j)
                            (1.d0 - dabs(s_y(mw,i,j))*dtdy) * wave_y_tilde(m,mw,i,j)
                    enddo
                    gp(m,i,j) = gp(m,i,j) + 0.5d0 * cqyy(m)
                    gm(m,i,j) = gm(m,i,j) + 0.5d0 * cqyy(m)
                enddo
            endif ! end if second-order 
    ! ##### solve for transverse waves and add to fp and fm
            if (method(3).ne.0) then ! if transverse propagation
                ! reconstruct bmdq and bpdq
                do m=1,meqn
                    bmdq(m) = s_y(1,i,j)*wave_y(m,1,i,j)
                    bpdq(m) = s_y(2,i,j)*wave_y(m,2,i,j)
                enddo
                if (method(2).gt.1 .and. method(3).eq.2) then
                    ! incorporate cqyy into bmdq and bpdq so that it is split also.
                    ! also reconstruct bmdq and bpdq
                    do m=1,meqn
                        bmdq(m) = bmdq(m) + cqyy(m)
                        bpdq(m) = bpdq(m) - cqyy(m)
                    enddo
                endif

                mu = 3
                mv = 2
                ! ##### solve for apbmdq and ambmdq
                a1 = (-bmdq(1) + zz*bmdq(mv)) / (2.d0*zz)
                a2 = (bmdq(1) + zz*bmdq(mv)) / (2.d0*zz)
                ambmdq(1) = cc * a1*zz
                ambmdq(mu) = 0.d0
                ambmdq(mv) = -cc * a1
                apbmdq(1) = cc * a2*zz
                apbmdq(mu) = 0.d0
                apbmdq(mv) = cc * a2

                if (mcapa > 0)  then
                    !todo
                    print * , "To be implemented"
                else
                    do m =1,meqn
                        fm(m,i,j-1) = fm(m,i,j-1) - 0.5d0*dtdy * ambmdq(m)
                        fp(m,i,j-1) = fp(m,i,j-1) - 0.5d0*dtdy * ambmdq(m)

                        fm(m,i+1,j-1) = fm(m,i+1,j-1) - 0.5d0*dtdy * apbmdq(m)
                        fp(m,i+1,j-1) = fp(m,i+1,j-1) - 0.5d0*dtdy * apbmdq(m)
                    enddo
                endif

                ! # solve for bpapdq and bmapdq
                a1 = (-bpdq(1) + zz*bpdq(mv)) / (2.d0*zz)
                a2 = (bpdq(1) + zz*bpdq(mv)) / (2.d0*zz)
                !        # The down-going flux difference bmasdq is the product  -c * wave
                ambpdq(1) = cc * a1*zz
                ambpdq(mu) = 0.d0
                ambpdq(mv) = -cc * a1
                !        # The up-going flux difference bpasdq is the product  c * wave
                apbpdq(1) = cc * a2*zz
                apbpdq(mu) = 0.d0
                apbpdq(mv) = cc * a2

                if (mcapa > 0)  then
                    !todo
                    print * , "To be implemented"
                else
                    do m =1,meqn
                        fm(m,i,j) = fm(m,i,j) - 0.5d0*dtdy * ambpdq(m)
                        fp(m,i,j) = fp(m,i,j) - 0.5d0*dtdy * ambpdq(m)

                        fm(m,i+1,j) = fm(m,i+1,j) - 0.5d0*dtdy * apbpdq(m)
                        fp(m,i+1,j) = fp(m,i+1,j) - 0.5d0*dtdy * apbpdq(m)
                    enddo
                endif
            endif ! end if transverse propagation
    ! ##### END solve for transverse waves and add to fp and fm
        enddo ! end my loop
    enddo ! end mx loop
end subroutine y_sweep_2nd_order


! For now, this kernel assumes 
! 1) limit == true
! 2) use 2nd order
! 3) include 2nd-order terms in transverse waves
! 4) van Leer limiter is used
! 5) use_fwaves == false
! 6) mcapa = 0
! We should write separate kernels for cases where 
! the options above are different
attributes(global) &
subroutine y_sweep_2nd_order_gpu(fm, fp, gm, gp, s_y, wave_y, &
    mbc, mx, my, dtdy, cc, zz)

    implicit none

    integer, value, intent(in) :: mbc, mx, my
    real(kind=8), intent(inout) :: fm(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: fp(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: gm(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: gp(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(in) ::    s_y(2-mbc:mx+mbc-1, 2-mbc:my + mbc, NWAVES)
    real(kind=8), intent(in) :: wave_y(2-mbc:mx+mbc-1, 2-mbc:my+mbc, NEQNS, NWAVES)
    real(kind=8), value, intent(in) :: dtdy
    real(kind=8), value, intent(in) :: cc, zz

    ! Local variables for the Riemann solver
    real(kind=8) :: cqyy(NEQNS)
    real(kind=8) :: bmdq(NEQNS), bpdq(NEQNS)
    real(kind=8) :: wave_y_tilde(NEQNS, NWAVES)
    real(kind=8) :: apbmdq(NEQNS), ambmdq(NEQNS), apbpdq(NEQNS), ambpdq(NEQNS)
    real(kind=8) :: a1, a2
    real(kind=8) :: dot, wnorm2, wlimitr, r
    integer :: i,j
    integer :: m, mw
    real(kind=8) :: atomic_result


    i = (blockIdx%x-1) * blockDim%x + threadIdx%x
    j = (blockIdx%y-1) * blockDim%y + threadIdx%y
    ! we shift i and j such that they are mapped to the loop:
    ! do i = 0, mx+1 ! mx loop
    !     do j = 1, my+1 ! my loop
    i = i + 0 - 1
    ! j is already corresponding to j = 1

    if (i > (mx+1) .or. j > (my+1) ) then
        return
    endif

    !     -----------------------------------------------------------
    !     # modify G fluxes for second order q_{yy} correction terms
    !     # and solve for transverse waves
    !     -----------------------------------------------------------
    do mw=1,NWAVES ! NWAVES loop
        dot = 0.d0
        wnorm2 = 0.d0
        do m=1,NEQNS
            wnorm2 = wnorm2 + wave_y(i,j,m,mw)**2
        enddo
        if (wnorm2.eq.0.d0) then
            wave_y_tilde(:,mw) =  wave_y(i,j,:,mw)
        else

            if (s_y(i,j,mw) .gt. 0.d0) then
                do m=1,NEQNS
                    dot = dot + wave_y(i,j,m,mw)*wave_y(i,j-1,m,mw)
                enddo
            else
                do m=1,NEQNS
                    dot = dot + wave_y(i,j,m,mw)*wave_y(i,j+1,m,mw)
                enddo
            endif

            r = dot / wnorm2

            !               ----------
            !               # van Leer
            !               ----------
            wlimitr = (r + dabs(r)) / (1.d0 + dabs(r))

            !
            !  # apply limiter to waves:
            !
            do m=1,NEQNS
                wave_y_tilde(:,mw) = wlimitr * wave_y(i,j,:,mw)
            enddo
        endif
    enddo ! end mwave loop

    cqyy = 0.d0
    do mw = 1, NWAVES
        cqyy(:) = cqyy(:) + &
            dabs(s_y(i,j,mw)) * (1.d0 - dabs(s_y(i,j,mw))*dtdy) * wave_y_tilde(:,mw)
    enddo


    gp(i,j,:) = gp(i,j,:) + 0.5d0 * cqyy(:)
    gm(i,j,:) = gm(i,j,:) + 0.5d0 * cqyy(:)
    ! reconstruct bmdq and bpdq
    bmdq(:) = s_y(i,j,1)*wave_y(i,j,:,1)
    bpdq(:) = s_y(i,j,2)*wave_y(i,j,:,2)

    bmdq(:) = bmdq(:) + cqyy(:)
    bpdq(:) = bpdq(:) - cqyy(:)

    ! ##### solve for apbmdq and ambmdq
    a1 = (-bmdq(1) + zz*bmdq(2)) / (2.d0*zz)
    a2 = (bmdq(1) + zz*bmdq(2)) / (2.d0*zz)
    ambmdq(1) = cc * a1*zz
    ambmdq(2) = -cc * a1
    ambmdq(3) = 0.d0
    apbmdq(1) = cc * a2*zz
    apbmdq(2) = cc * a2
    apbmdq(3) = 0.d0

    do m =1,NEQNS
        atomic_result = atomicadd(fm(i,j-1  ,m), - 0.5d0*dtdy * ambmdq(m))
        atomic_result = atomicadd(fp(i,j-1  ,m), - 0.5d0*dtdy * ambmdq(m))
        atomic_result = atomicadd(fm(i+1,j-1,m), - 0.5d0*dtdy * apbmdq(m))
        atomic_result = atomicadd(fp(i+1,j-1,m), - 0.5d0*dtdy * apbmdq(m))
    enddo

    ! # solve for apbpdq and ambpdq
    a1 = (-bpdq(1) + zz*bpdq(2)) / (2.d0*zz)
    a2 = (bpdq(1) + zz*bpdq(2)) / (2.d0*zz)
    !        # The down-going flux difference bmasdq is the product  -c * wave
    ambpdq(1) = cc * a1*zz
    ambpdq(2) = -cc * a1
    ambpdq(3) = 0.d0
    !        # The up-going flux difference bpasdq is the product  c * wave
    apbpdq(1) = cc * a2*zz
    apbpdq(2) = cc * a2
    apbpdq(3) = 0.d0

    do m =1,NEQNS
        atomic_result = atomicadd(fm(i,j  ,m), - 0.5d0*dtdy * ambpdq(m))
        atomic_result = atomicadd(fp(i,j  ,m), - 0.5d0*dtdy * ambpdq(m))
        atomic_result = atomicadd(fm(i+1,j,m), - 0.5d0*dtdy * apbpdq(m))
        atomic_result = atomicadd(fp(i+1,j,m), - 0.5d0*dtdy * apbpdq(m))
    enddo
    return
end subroutine y_sweep_2nd_order_gpu

end module sweep_module
