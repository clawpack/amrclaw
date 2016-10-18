subroutine step1(maxm,meqn,maux,mbc,mx,qold,aux,dx,dt,cflgrid,fm,fp,rp1)
!
!     clawpack routine ...  modified for AMRCLAW
!
!     Take one time step, updating q.
!     On entry, qold gives
!        initial data for this step
!        and is unchanged in this version.
!    
!     fm, fp are fluxes to left and right of single cell edge
!     See the flux1 documentation for more information.
!
!     Converted to f90 2012-1-04 (KTM)
!     Modified from step2 to work with 1d 2016-6-6 (BND)
    
    use amr_module

    implicit none
    
    external rp1
    
    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx
    real(kind=8), intent(in) :: dx,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: qold(meqn, 1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc)
    real(kind=8), intent(inout) :: fp(meqn,1-mbc:mx+mbc)
    
    ! Local storage for flux accumulation
    real(kind=8) :: faddm(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: faddp(meqn,1-mbc:maxm+mbc)
    
    ! Scratch storage for Sweeps and Riemann problems
    real(kind=8) ::  q1d(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: aux1(maux,1-mbc:maxm+mbc)
    real(kind=8) :: dtdx1d(1-mbc:maxm+mbc)
    real(kind=8) :: dtdy1d(1-mbc:maxm+mbc)
    
    real(kind=8) ::  wave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8) ::     s(mwaves, 1-mbc:maxm + mbc)
    real(kind=8) ::  amdq(meqn,1-mbc:maxm + mbc)
    real(kind=8) ::  apdq(meqn,1-mbc:maxm + mbc)
    real(kind=8) ::  cqxx(meqn,1-mbc:maxm + mbc)
    real(kind=8) :: bmadq(meqn,1-mbc:maxm + mbc)
    real(kind=8) :: bpadq(meqn,1-mbc:maxm + mbc)
    
    ! Looping scalar storage
    integer :: i,j,thread_num
    real(kind=8) :: dtdx,cfl1d
    
    ! Common block storage
    integer :: icom
    real(kind=8) :: dtcom,dxcom,tcom
    common /comxt/ dtcom,dxcom,tcom,icom
    
    ! Store mesh parameters in common block
    dxcom = dx
    dtcom = dt
    
    cflgrid = 0.d0
    dtdx = dt/dx
    
    fm = 0.d0
    fp = 0.d0
    
    ! ============================================================================
    ! Perform X-Sweep
        ! Copy old q into 1d slice
        q1d(:,1-mbc:mx+mbc) = qold(:,1-mbc:mx+mbc)
        
        ! Set dtdx slice if a capacity array exists
        if (mcapa > 0)  then
            dtdx1d(1-mbc:mx+mbc) = dtdx / aux(mcapa,1-mbc:mx+mbc)
        else
            dtdx1d = dtdx
        endif
        
        ! Copy aux array into slices
        if (maux > 0) then
            aux1(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc)
        endif

        ! Compute modification fadd to fluxes along this slice:
        call flux1(maxm,meqn,maux,mbc,mx,q1d,dtdx1d,aux1, &
                   faddm,faddp,cfl1d,wave,s, &
                   amdq,apdq,cqxx,bmadq,bpadq,rp1)
                   
        cflgrid = max(cflgrid,cfl1d)

        ! Update fluxes
        fm(:,1:mx+1) = fm(:,1:mx+1) + faddm(:,1:mx+1)
        fp(:,1:mx+1) = fp(:,1:mx+1) + faddp(:,1:mx+1)

end subroutine step1
