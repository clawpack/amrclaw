!
! -------------------------------------------------------------
!
subroutine step2x(maxm,meqn,maux,mbc,mx,my,qold,aux,dx,dt,cflgrid,fm,fp,rpn2)
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
    real(kind=8), intent(in) :: dx,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

   
    ! Local storage for flux accumulation
    real(kind=8) :: faddm(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: faddp(meqn,1-mbc:maxm+mbc)
 

    ! Scratch storage for Sweeps and Riemann problems
    real(kind=8) ::  q1d(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: aux2(maux,1-mbc:maxm+mbc)
    real(kind=8) :: dtdx1d(1-mbc:maxm+mbc)

    
    real(kind=8) ::  wave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8) ::     s(mwaves, 1-mbc:maxm + mbc)
    real(kind=8) ::  cqxx(meqn,1-mbc:maxm + mbc)
   
    
    ! Looping scalar storage
    integer :: j
    real(kind=8) :: dtdx,cfl1d
    
    
    cflgrid = 0.d0
    dtdx = dt/dx
   
    fm = 0.d0
    fp = 0.d0
    
    ! ============================================================================
    ! Perform X-Sweeps (for Godunov split method, called from stepgrid_dimSplit)
    ! loop indices assume x sweep goes first
    do j = 1-mbc,my+mbc   

        ! Copy old q into 1d slice
        q1d(:,1-mbc:mx+mbc) = qold(:,1-mbc:mx+mbc,j)
        
        ! Set dtdx slice if a capacity array exists
        if (mcapa > 0)  then
            dtdx1d(1-mbc:mx+mbc) = dtdx / aux(mcapa,1-mbc:mx+mbc,j)
        else
            dtdx1d = dtdx
        endif
        
        ! Copy aux array into slices
        if (maux > 0) then

            aux2(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j  )

        endif
        
        ! Compute modifications fadd and gadd to fluxes along this slice:
!!!     ::: dont need aux1 and aux3 for dimensionally split method
!!!     ::: but flux2 needs 3 aux arrays, so passing in aux2 3 times
        call flux2_dimSplit(1,maxm,meqn,maux,mbc,mx,q1d,dtdx1d,aux2, &
                   faddm,faddp,cfl1d,wave,s,cqxx,rpn2)       
                   
        cflgrid = max(cflgrid,cfl1d)

        ! Update fluxes
        fm(:,1:mx+1,j) = fm(:,1:mx+1,j) + faddm(:,1:mx+1)
        fp(:,1:mx+1,j) = fp(:,1:mx+1,j) + faddp(:,1:mx+1)
         
    enddo


end subroutine step2x
