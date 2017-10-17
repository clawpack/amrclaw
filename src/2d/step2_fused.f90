!> Compute all fluxes at cell edges 
!! \param qold[in] solution array for computing fluxes. It is not changed in this subroutine
!! \param fm[out] fluxes on the left side of each vertical edge
!! \param fp[out] fluxes on the right side of each vertical edge
!! \param gm[out] fluxes on the lower side of each horizontal edge
!! \param gp[out] fluxes on the upper side of each horizontal edge
subroutine step2_fused(maxm,meqn,maux,mbc,mx,my,qold,aux,dx,dy,dt,cflgrid,fm,fp,gm,gp,rpn2,rpt2)
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

    implicit none
    
    external rpn2, rpt2
    
    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my
    real(kind=8), intent(in) :: dx,dy,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    
    ! Local storage for flux accumulation
    real(kind=8) :: faddm(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: faddp(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: gaddm(meqn,1-mbc:maxm+mbc,2)
    real(kind=8) :: gaddp(meqn,1-mbc:maxm+mbc,2)
    
    ! Scratch storage for Sweeps and Riemann problems
    real(kind=8) :: aux1(maux,1-mbc:maxm+mbc)
    real(kind=8) :: aux2(maux,1-mbc:maxm+mbc)
    real(kind=8) :: aux3(maux,1-mbc:maxm+mbc)
    
    real(kind=8) :: bmadq(meqn,1-mbc:maxm + mbc)
    real(kind=8) :: bpadq(meqn,1-mbc:maxm + mbc)
    
    ! Looping scalar storage
    integer :: i,j,thread_num
    real(kind=8) :: dtdx,dtdy

    ! Local variables for the Riemann solver
    real(kind=8) :: delta1, delta2, a1, a2
    integer :: m, mw, mu, mv
    real(kind=8) :: rho, bulk, cc, zz
    real(kind=8) :: s_x(mwaves, 1-mbc:mx + mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: s_y(mwaves, 2-mbc:mx+mbc-1, 1-mbc:my + mbc)
    real(kind=8) :: wave_x(meqn, mwaves, 1-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: wave_y(meqn, mwaves, 2-mbc:mx+mbc-1, 1-mbc:my+mbc)
    real(kind=8) :: dtdxr, dtdxl, dtdyr, dtdyl
    real(kind=8) :: amdq, apdq
    real(kind=8) :: bmdq(meqn,1-mbc:my + mbc)
    real(kind=8) :: bpdq(meqn,1-mbc:my + mbc)
    ! For 2nd order corrections
    real(kind=8) :: wave_x_tmp(meqn, mwaves, 1-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: wave_y_tmp(meqn, mwaves, 2-mbc:mx+mbc-1, 1-mbc:my+mbc)
    real(kind=8) :: dot, wnorm2, wlimitr, dtdxave, dtdyave, abs_sign, c, r
    real(kind=8) :: cqxx(meqn,1-mbc:mx + mbc)
    real(kind=8) :: cqyy(meqn,1-mbc:my + mbc)
    logical limit

    common /cparam/ rho,bulk,cc,zz

    
    ! Common block storage
    ! integer :: icom,jcom
    ! real(kind=8) :: dtcom,dxcom,dycom,tcom
    ! common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
    
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

    limit = .false.
    do mw=1,mwaves
        if (mthlim(mw) .gt. 0) limit = .true.
    enddo
    
    ! ============================================================================
    ! Perform X-Sweeps
    do j = 0,my+1
        do i = 2-mbc, mx+mbc
            ! solve Riemann problem between cell (i-1,j) and (i,j)
            mu = 2
            mv = 3
            delta1 = qold(1,i,j) - qold(1,i-1,j)
            delta2 = qold(mu,i,j) - qold(mu,i-1,j)
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
                amdq = s_x(1,i,j)*wave_x(m,1,i,j)
                apdq = s_x(2,i,j)*wave_x(m,2,i,j)
                fm(m,i,j) = fm(m,i,j) + amdq
                fp(m,i,j) = fp(m,i,j) - apdq
            enddo
            if (mcapa > 0)  then
                dtdxl = dtdx / aux(mcapa,i-1,j)
                dtdxr = dtdx / aux(mcapa,i,j)
            else
                dtdxl = dtdx
                dtdxr = dtdx
            endif
            do mw=1,mwaves
                cflgrid = dmax1(cflgrid, dtdxr*s_x(mw,i,j),-dtdxl*s_x(mw,i,j))
            enddo
        enddo
    enddo
!     -----------------------------------------------------------
!     # modify F fluxes for second order q_{xx} correction terms:
!     -----------------------------------------------------------
    if (method(2).ne.1) then ! if second-order
        wave_x_tmp = wave_x
        do j = 0,my+1 ! my loop
            do i = 1, mx+1 ! mx loop
                ! # apply limiter to waves:
                if (limit) then ! limiter if
                    do mw=1,mwaves ! mwaves loop
                        if (mthlim(mw) .eq. 0) cycle
                        dot = 0.d0
                        wnorm2 = 0.d0
                        do m=1,meqn
                            wnorm2 = wnorm2 + wave_x_tmp(m,mw,i,j)**2
                        enddo
                        if (wnorm2.eq.0.d0) cycle

                        if (s_x(mw,i,j) .gt. 0.d0) then
                            do m=1,meqn
                                dot = dot + wave_x_tmp(m,mw,i,j)*wave_x_tmp(m,mw,i-1,j)
                            enddo
                        else
                            do m=1,meqn
                                dot = dot + wave_x_tmp(m,mw,i,j)*wave_x_tmp(m,mw,i+1,j)
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
                            wave_x(m,mw,i,j) = wlimitr * wave_x(m,mw,i,j)
                        enddo
                    enddo ! end mwave loop
                endif ! end limiter if

                if (mcapa > 0)  then
                    dtdxl = dtdx / aux(mcapa,i-1,j)
                    dtdxr = dtdx / aux(mcapa,i,j)
                else
                    dtdxl = dtdx
                    dtdxr = dtdx
                endif
                dtdxave = 0.5d0 * (dtdxl + dtdxr)
                ! second order corrections:
                do m=1,meqn
                    cqxx(m,i) = 0.d0
                    do mw=1,mwaves
                        if (use_fwaves) then
                            abs_sign = dsign(1.d0,s_x(mw,i,j))
                        else
                            abs_sign = dabs(s_x(mw,i,j))
                        endif

                        cqxx(m,i) = cqxx(m,i) + abs_sign * &
                            (1.d0 - dabs(s_x(mw,i,j))*dtdxave) * wave_x(m,mw,i,j)
                    enddo
                    fp(m,i,j) = fp(m,i,j) + 0.5d0 * cqxx(m,i)
                    fm(m,i,j) = fm(m,i,j) + 0.5d0 * cqxx(m,i)
                enddo
            enddo ! end mx loop
        enddo ! end my loop
    endif ! end if second-order 
!     -----------------------------------------------------------
!     # END modify F fluxes for second order q_{xx} correction terms:
!     -----------------------------------------------------------

        ! if (method(3).ne.0) then !# has transverse propagation
        !     if (method(2).gt.1 .and. method(3).eq.2) then
        !         ! # incorporate cqxx into amdq and apdq so that it is split also.
        !         do i = 1, mx+1
        !             do m=1,meqn
        !                 amdq(m,i) = amdq(m,i) + cqxx(m,i)
        !                 apdq(m,i) = apdq(m,i) - cqxx(m,i)
        !             enddo
        !         enddo
        !     endif
        ! endif

!     -----------------------------------------------------------
!     # modify G fluxes for transverse propagation
!     -----------------------------------------------------------
        ! # split the left-going flux difference into down-going and up-going:
        ! # split the right-going flux difference into down-going and up-going:
!     -----------------------------------------------------------
!     # END modify G fluxes for transverse propagation
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
                bmdq(m,j) = s_y(1,i,j)*wave_y(m,1,i,j)
                bpdq(m,j) = s_y(2,i,j)*wave_y(m,2,i,j)
                gm(m,i,j) = gm(m,i,j) + bmdq(m,j)
                gp(m,i,j) = gp(m,i,j) - bpdq(m,j)
            enddo
            if (mcapa > 0)  then
                dtdyl = dtdy / aux(mcapa,i-1,j)
                dtdyr = dtdy / aux(mcapa,i,j)
            else
                dtdyl = dtdy
                dtdyr = dtdy
            endif
            do mw=1,mwaves
                cflgrid = dmax1(cflgrid, dtdyr*s_y(mw,i,j),-dtdyl*s_y(mw,i,j))
            enddo
        enddo
    enddo
!     -----------------------------------------------------------
!     # modify G fluxes for second order q_{yy} correction terms:
!     -----------------------------------------------------------
    if (method(2).ne.1) then ! if second-order
        wave_y_tmp = wave_y
        do i = 0, mx+1 ! mx loop
            do j = 1, my+1 ! my loop
                ! # apply limiter to waves:
                if (limit) then ! limiter if
                    do mw=1,mwaves ! mwaves loop
                        if (mthlim(mw) .eq. 0) cycle
                        dot = 0.d0
                        wnorm2 = 0.d0
                        do m=1,meqn
                            wnorm2 = wnorm2 + wave_y_tmp(m,mw,i,j)**2
                        enddo
                        if (wnorm2.eq.0.d0) cycle

                        if (s_y(mw,i,j) .gt. 0.d0) then
                            do m=1,meqn
                                dot = dot + wave_y_tmp(m,mw,i,j)*wave_y_tmp(m,mw,i,j-1)
                            enddo
                        else
                            do m=1,meqn
                                dot = dot + wave_y_tmp(m,mw,i,j)*wave_y_tmp(m,mw,i,j+1)
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
                            wave_y(m,mw,i,j) = wlimitr * wave_y(m,mw,i,j)
                        enddo
                    enddo ! end mwave loop
                endif ! end limiter if
                if (mcapa > 0)  then
                    dtdyl = dtdy / aux(mcapa,i,j-1)
                    dtdyr = dtdy / aux(mcapa,i,j)
                else
                    dtdyl = dtdy
                    dtdyr = dtdy
                endif
                dtdyave = 0.5d0 * (dtdyl + dtdyr)
                ! second order corrections:
                do m=1,meqn
                    cqyy(m,j) = 0.d0
                    do mw=1,mwaves
                        if (use_fwaves) then
                            abs_sign = dsign(1.d0,s_y(mw,i,j))
                        else
                            abs_sign = dabs(s_y(mw,i,j))
                        endif

                        cqyy(m,j) = cqyy(m,j) + abs_sign * &
                            (1.d0 - dabs(s_y(mw,i,j))*dtdyave) * wave_y(m,mw,i,j)
                    enddo
                    gp(m,i,j) = gp(m,i,j) + 0.5d0 * cqyy(m,j)
                    gm(m,i,j) = gm(m,i,j) + 0.5d0 * cqyy(m,j)
                enddo
            enddo ! end my loop
        enddo ! end mx loop
    endif ! end if second-order 
!     -----------------------------------------------------------
!     # END modify G fluxes for second order q_{yy} correction terms:
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
