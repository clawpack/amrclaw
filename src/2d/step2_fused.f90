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
    real(kind=8) ::  q1d(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: aux1(maux,1-mbc:maxm+mbc)
    real(kind=8) :: aux2(maux,1-mbc:maxm+mbc)
    real(kind=8) :: aux3(maux,1-mbc:maxm+mbc)
    real(kind=8) :: dtdx1d(1-mbc:maxm+mbc)
    real(kind=8) :: dtdy1d(1-mbc:maxm+mbc)
    
    real(kind=8) ::  wave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8) ::     s(mwaves, 1-mbc:maxm + mbc)
    real(kind=8) ::  cqxx(meqn,1-mbc:maxm + mbc)
    real(kind=8) :: bmadq(meqn,1-mbc:maxm + mbc)
    real(kind=8) :: bpadq(meqn,1-mbc:maxm + mbc)
    
    ! Looping scalar storage
    integer :: i,j,thread_num
    real(kind=8) :: dtdx,dtdy,cfl1d

    ! Local variables for the Riemann solver
    real(kind=8) :: delta1, delta2, a1, a2
    integer :: m, mw, mu, mv
    real(kind=8) :: rho, bulk, cc, zz
    real(kind=8) ::  amdq(meqn,1-mbc:mx + mbc)
    real(kind=8) ::  apdq(meqn,1-mbc:mx + mbc)
    real(kind=8) ::  bmdq(meqn,1-mbc:my + mbc)
    real(kind=8) ::  bpdq(meqn,1-mbc:my + mbc)

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
    
    ! ============================================================================
    ! Perform X-Sweeps
    do j = 0,my+1

        ! Set dtdx slice if a capacity array exists
        if (mcapa > 0)  then
            dtdx1d(1-mbc:mx+mbc) = dtdx / aux(mcapa,1-mbc:mx+mbc,j)
        else
            dtdx1d = dtdx
        endif

        ! do i = 2-mbc, mx+mbc
        do i = 1, mx+1
            ! solve Riemann problem between cell (i-1,j) and (i,j)
            mu = 2
            mv = 3
            delta1 = qold(1,i,j) - qold(1,i-1,j)
            delta2 = qold(mu,i,j) - qold(mu,i-1,j)
            a1 = (-delta1 + zz*delta2) / (2.d0*zz)
            a2 = (delta1 + zz*delta2) / (2.d0*zz)
            !        # Compute the waves.
            wave(1,1,i) = -a1*zz
            wave(mu,1,i) = a1
            wave(mv,1,i) = 0.d0
            s(1,i) = -cc

            wave(1,2,i) = a2*zz
            wave(mu,2,i) = a2
            wave(mv,2,i) = 0.d0
            s(2,i) = cc
            do m = 1,meqn
                amdq(m,i) = s(1,i)*wave(m,1,i)
                apdq(m,i) = s(2,i)*wave(m,2,i)
                fm(m,i,j) = fm(m,i,j) + amdq(m,i)
                fp(m,i,j) = fp(m,i,j) - apdq(m,i)
            enddo
            do mw=1,mwaves
                cflgrid = dmax1(cflgrid, dtdx1d(i)*s(mw,i),-dtdx1d(i-1)*s(mw,i))
            enddo
        enddo
    enddo


    ! ============================================================================
    !  y-sweeps    
    !

    do i = 0,mx+1

        ! Set dtdx slice if a capacity array exists
        if (mcapa > 0) then
            dtdy1d(1-mbc:my+mbc) = dtdy / aux(mcapa,i,1-mbc:my+mbc)
        else
            dtdy1d = dtdy
        endif

        ! do j = 2-mbc, my+mbc
        do j = 1, my+1
            ! solve Riemann problem between cell (i-1,j) and (i,j)
            mu = 3
            mv = 2
            delta1 = qold(1,i,j) - qold(1,i,j-1)
            delta2 = qold(mu,i,j) - qold(mu,i,j-1)
            a1 = (-delta1 + zz*delta2) / (2.d0*zz)
            a2 = (delta1 + zz*delta2) / (2.d0*zz)
            !        # Compute the waves.
            wave(1,1,j) = -a1*zz
            wave(mu,1,j) = a1
            wave(mv,1,j) = 0.d0
            s(1,j) = -cc

            wave(1,2,j) = a2*zz
            wave(mu,2,j) = a2
            wave(mv,2,j) = 0.d0
            s(2,j) = cc
            do m = 1,meqn
                bmdq(m,j) = s(1,j)*wave(m,1,j)
                bpdq(m,j) = s(2,j)*wave(m,2,j)
                gm(m,i,j) = gm(m,i,j) + bmdq(m,j)
                gp(m,i,j) = gp(m,i,j) - bpdq(m,j)
            enddo
            do mw=1,mwaves
                cflgrid = dmax1(cflgrid, dtdy1d(j)*s(mw,j),-dtdy1d(j-1)*s(mw,j))
            enddo
        enddo
    enddo


end subroutine step2_fused
