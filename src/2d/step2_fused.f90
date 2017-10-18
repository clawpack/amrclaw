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

    real(kind=8) :: s_x(mwaves, 1-mbc:mx + mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: s_y(mwaves, 2-mbc:mx+mbc-1, 1-mbc:my + mbc)
    real(kind=8) :: wave_x(meqn, mwaves, 1-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: wave_y(meqn, mwaves, 2-mbc:mx+mbc-1, 1-mbc:my+mbc)


    
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
    
!     -----------------------------------------------------------
!     # modify F fluxes for second order q_{xx} correction terms
!     # and solve for transverse waves
!     -----------------------------------------------------------

    call x_sweep_2nd_order(fm, fp, gm, gp, s_x, wave_x, meqn, mwaves, mbc, mx, my, dtdx)


    call y_sweep_1st_order(qold, gm, gp, s_y, wave_y, meqn, mwaves, mbc, mx, my, dtdy, cflgrid)

!     -----------------------------------------------------------
!     # modify G fluxes for second order q_{yy} correction terms:
!     # and transverse waves
!     -----------------------------------------------------------

    call y_sweep_2nd_order(fm, fp, gm, gp, s_y, wave_y, meqn, mwaves, mbc, mx, my, dtdy)

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
