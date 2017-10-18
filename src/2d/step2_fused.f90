!> Compute all fluxes at cell edges 
!! \param qold[in] solution array for computing fluxes. It is not changed in this subroutine
!! \param fm[out] fluxes on the left side of each vertical edge
!! \param fp[out] fluxes on the right side of each vertical edge
!! \param gm[out] fluxes on the lower side of each horizontal edge
!! \param gp[out] fluxes on the upper side of each horizontal edge
subroutine step2_fused(maxm,meqn,maux,mbc,mx,my,qold,dx,dy,dt,cflgrid,fm,fp,gm,gp,rpn2,rpt2)
!
!     clawpack routine ...  modified for AMRCLAW
!
!     Take one time step, updating q.
!     On entry, qold gives
!        initial data for this step
!        and is unchanged in this version.
!    
!     fm, fp are fluxes to left and right of single cell edge
!     See the flux2 documentation for more information.
!
!     Converted to f90 2012-1-04 (KTM)
!
    
    use amr_module
    use parallel_advanc_module, only: dtcom, dxcom, dycom, icom, jcom
    use sweep_module, only: x_sweep_1st_order, x_sweep_2nd_order, y_sweep_1st_order, y_sweep_2nd_order 

    implicit none
    
    external rpn2, rpt2
    
    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my
    real(kind=8), intent(in) :: dx,dy,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    ! real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

#ifdef CUDA
    attributes(device) :: qold, fm, fp, gm, gp
#endif
    
    ! Looping scalar storage
    integer :: i,j
    real(kind=8) :: dtdx,dtdy

    ! Local variables for the Riemann solver
    real(kind=8) :: delta1, delta2, a1, a2
    integer :: m, mw, mu, mv
    real(kind=8) :: rho, bulk, cc, zz
    real(kind=8) :: s_x(mwaves, 1-mbc:mx + mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: s_y(mwaves, 2-mbc:mx+mbc-1, 1-mbc:my + mbc)
    real(kind=8) :: wave_x(meqn, mwaves, 1-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: wave_y(meqn, mwaves, 2-mbc:mx+mbc-1, 1-mbc:my+mbc)
    real(kind=8) :: amdq(meqn), apdq(meqn)
    real(kind=8) :: bmdq(meqn), bpdq(meqn)

    ! For 2nd order corrections
    real(kind=8) :: wave_x_tilde(meqn, mwaves, 1-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: wave_y_tilde(meqn, mwaves, 2-mbc:mx+mbc-1, 1-mbc:my+mbc)

    ! real(kind=8) :: dtdxave, dtdyave
    real(kind=8) :: dot, wnorm2, wlimitr, abs_sign, c, r
    real(kind=8) :: cqxx(meqn)
    real(kind=8) :: cqyy(meqn)

    ! For transverse waves
    real(kind=8) :: bpamdq(meqn), bmamdq(meqn), bpapdq(meqn), bmapdq(meqn)
    real(kind=8) :: apbmdq(meqn), ambmdq(meqn), apbpdq(meqn), ambpdq(meqn)
    logical limit

    common /cparam/ rho,bulk,cc,zz

    
    ! Store mesh parameters in common block
    dxcom = dx
    dycom = dy
    dtcom = dt
    
    cflgrid = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy

    fm = 0.d0
    fp = 0.d0
    gm = 0.d0
    gp = 0.d0

    call x_sweep_1st_order(qold, fm, fp, s_x, wave_x, meqn, mwaves, mbc, mx, my, dtdx, cflgrid)
    ! call x_sweep_2nd_order(qold, fm, fp, gm, gp, s_x, wave_x)
    ! call y_sweep_1st_order(qold, gm, gp, s_x, wave_x)
    ! call y_sweep_2nd_order(qold, fm, fp, gm, gp, s_x, wave_x)

    limit = .false.
    do mw=1,mwaves
        if (mthlim(mw) .gt. 0) limit = .true.
    enddo
    
    ! ! ============================================================================
    ! ! Perform X-Sweeps
    ! do j = 0,my+1
    !     do i = 2-mbc, mx+mbc
    !         ! solve Riemann problem between cell (i-1,j) and (i,j)
    !         mu = 2
    !         mv = 3
    !         delta1 = qold(1,i,j) - qold(1,i-1,j)
    !         delta2 = qold(mu,i,j) - qold(mu,i-1,j)
    !         a1 = (-delta1 + zz*delta2) / (2.d0*zz)
    !         a2 = (delta1 + zz*delta2) / (2.d0*zz)
    !         !        # Compute the waves.
    !         wave_x(1,1,i,j) = -a1*zz
    !         wave_x(mu,1,i,j) = a1
    !         wave_x(mv,1,i,j) = 0.d0
    !         s_x(1,i,j) = -cc

    !         wave_x(1,2,i,j) = a2*zz
    !         wave_x(mu,2,i,j) = a2
    !         wave_x(mv,2,i,j) = 0.d0
    !         s_x(2,i,j) = cc
    !         do m = 1,meqn
    !             amdq(m) = s_x(1,i,j)*wave_x(m,1,i,j)
    !             apdq(m) = s_x(2,i,j)*wave_x(m,2,i,j)
    !             if (i >= 1 .and. i<=(mx+1)) then
    !                 fm(m,i,j) = fm(m,i,j) + amdq(m)
    !                 fp(m,i,j) = fp(m,i,j) - apdq(m)
    !             endif
    !         enddo
    !         ! if (mcapa > 0)  then
    !         !     dtdxl = dtdx / aux(mcapa,i-1,j)
    !         !     dtdxr = dtdx / aux(mcapa,i,j)
    !         ! else
    !         !     dtdxl = dtdx
    !         !     dtdxr = dtdx
    !         ! endif
    !         do mw=1,mwaves
    !             if (i >= 1 .and. i<=(mx+1)) then
    !                 ! cflgrid = dmax1(cflgrid, dtdxr*s_x(mw,i,j),-dtdxl*s_x(mw,i,j))
    !                 cflgrid = dmax1(cflgrid, dtdx*s_x(mw,i,j),-dtdx*s_x(mw,i,j))
    !             endif
    !         enddo
    !     enddo
    ! enddo
!     -----------------------------------------------------------
!     # modify F fluxes for second order q_{xx} correction terms
!     # and solve for transverse waves
!     -----------------------------------------------------------
    if (method(2).ne.1) then ! if second-order
        wave_x_tilde = wave_x
    endif
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
!     -----------------------------------------------------------
!     # END modify F fluxes for second order q_{xx} correction terms:
!     # and solve for transverse waves
!     -----------------------------------------------------------


    ! ============================================================================
    !  y-sweeps    
    do i = 0,mx+1
        do j = 2-mbc, my+mbc
            ! solve Riemann problem between cell (i,j-1) and (i,j)
            mu = 3
            mv = 2
            delta1 = qold(1,i,j) - qold(1,i,j-1)
            delta2 = qold(mu,i,j) - qold(mu,i,j-1)
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
!     -----------------------------------------------------------
!     # modify G fluxes for second order q_{yy} correction terms:
!     # and transverse waves
!     -----------------------------------------------------------
    if (method(2).ne.1) then ! if second-order
        wave_y_tilde = wave_y
    endif
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
!     -----------------------------------------------------------
!     # END modify G fluxes for second order q_{yy} correction terms:
!     # and transverse waves
!     -----------------------------------------------------------


end subroutine step2_fused

! Algorithm1 (don't use shared memory):
!
! # x-sweep 
!
! kernel1:
! * READ q from DRAM
! Compute amdq and apdq:
!   store amdq, apdq, s and wave in register
! * READ fp and fm from DRAM
! add amdq and apdq to fp and fm
! * WRITE s, wave, fm, fp to DRAM
!
! kernel2:
! * READ s, wave, fm, fp from DRAM
! Use s and wave to compute cqxx
! Add cqxx to fm and fp
! Compute new amdq and apdq from s, wave and cqxx
! * READ gm and gp from DRAM
! Compute bpapdq, bmapdq, bpamdq, bmamdq and add them to gm and gp
! * WRITE gm and gp to DRAM
! 
! 
! # y-sweep 
!
! kernel1:
! * READ q from DRAM
! Compute bmdq and bpdq:
!   store bmdq, bpdq, s and wave in register
! * READ gp and gm from DRAM
! add bmdq and bpdq to gp and gm
! * WRITE s, wave, gm, gp to DRAM
!
! kernel2:
! * READ s, wave, gm, gp from DRAM
! Use s and wave to compute cqyy
! Add cqyy to gm and gp
! Compute new bmdq and bpdq from s, wave and cqyy
! * READ fm and fp from DRAM
! Compute apbpdq, ambpdq, apbmdq, ambmdq and add them to fm and fp
! * WRITE fm and fp to DRAM
! 


! Algorithm2 (use shared memory):
!
! # x-sweep 
!
! kernel:
! * READ q from DRAM
! Compute amdq and apdq:
!   store amdq, apdq, s and wave in register
! * READ fp and fm from DRAM
! add amdq and apdq to fp and fm
! * store s, wave in shared memory
! synchronize all threads in a block
! Use s and wave to compute cqxx
! Add cqxx to fm and fp
! Compute new amdq and apdq from s, wave and cqxx
! * READ gm and gp from DRAM
! Compute bpapdq, bmapdq, bpamdq, bmamdq and add them to gm and gp
! * WRITE fm, fp, gm and gp to DRAM
! 
! 
! # y-sweep 
!
! kernel1:
! * READ q from DRAM
! Compute bmdq and bpdq:
!   store bmdq, bpdq, s and wave in register
! * READ gp and gm from DRAM
! add bmdq and bpdq to gp and gm
! * store s, wave in shared memory
! synchronize all threads in a block
! Use s and wave to compute cqyy
! Add cqyy to gm and gp
! Compute new bmdq and bpdq from s, wave and cqyy
! * READ fm and fp from DRAM
! Compute apbpdq, ambpdq, apbmdq, ambmdq and add them to fm and fp
! * WRITE fm, fp, gm and gp to DRAM
