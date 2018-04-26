!
! -------------------------------------------------------------
!
subroutine step2y(maxm,meqn,maux,mbc,mx,my,qold,aux,dy,dt,cflgrid,gm,gp,rpn2)
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

    implicit none
    
    external rpn2
    
    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my
    real(kind=8), intent(in) :: dy,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    
    ! Local storage for flux accumulation (these help compute normal g flux, based on step2 usage)
    real(kind=8) :: faddm(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: faddp(meqn,1-mbc:maxm+mbc)

    
    ! Scratch storage for Sweeps and Riemann problems
    real(kind=8) ::  q1d(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: aux2(maux,1-mbc:maxm+mbc)
    real(kind=8) :: dtdy1d(1-mbc:maxm+mbc)
    
    real(kind=8) ::  wave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8) ::     s(mwaves, 1-mbc:maxm + mbc)
    real(kind=8) ::  cqxx(meqn,1-mbc:maxm + mbc)
    
    ! Looping scalar storage
    integer :: i
    real(kind=8) :: dtdy,cfl1d
    
    
    cflgrid = 0.d0
    dtdy = dt/dy
    
    gm = 0.d0
    gp = 0.d0


    ! ============================================================================
    !  y-sweeps   (for Godunov split method, called from stepgrid_dimSplit)
    !  loop indices assume that x sweep has gone first, this is second sweep.
    do i = 1,mx
        
        ! Copy data along a slice into 1d arrays:
        q1d(:,1-mbc:my+mbc) = qold(:,i,1-mbc:my+mbc)

        ! Set dt/dy ratio in slice
        if (mcapa > 0) then
            dtdy1d(1-mbc:my+mbc) = dtdy / aux(mcapa,i,1-mbc:my+mbc)
        else
            dtdy1d = dtdy
        endif

        ! Copy aux slices
        if (maux .gt. 0)  then

            aux2(:,1-mbc:my+mbc) = aux(:,i,1-mbc:my+mbc)

        endif
        
        ! Compute modifications fadd and gadd to fluxes along this slice
!!!     ::: dont need aux1 and aux3 for dimensionally split method
!!!     ::: but flux2 needs 3 aux arrays, so passing in aux2 3 times
        call flux2_dimSplit(2,maxm,meqn,maux,mbc,my,q1d,dtdy1d,aux2, &
                   faddm,faddp,cfl1d,wave,s,cqxx,rpn2)

        cflgrid = max(cflgrid,cfl1d)

        ! Update fluxes
        gm(:,i,1:my+1) = gm(:,i,1:my+1) + faddm(:,1:my+1)
        gp(:,i,1:my+1) = gp(:,i,1:my+1) + faddp(:,1:my+1)

    enddo


end subroutine step2y
