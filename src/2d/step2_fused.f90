!> Compute all fluxes at cell edges 
!! \param q[in] solution array for computing fluxes. It is not changed in this subroutine
!! \param fm[out] fluxes on the left side of each vertical edge
!! \param fp[out] fluxes on the right side of each vertical edge
!! \param gm[out] fluxes on the lower side of each horizontal edge
!! \param gp[out] fluxes on the upper side of each horizontal edge
subroutine step2_fused(maxm,meqn,maux,mbc,mx,my,q,dx,dy,dt,cflgrid,fm,fp,gm,gp,rpn2,rpt2)
!
!     clawpack routine ...  modified for AMRCLAW
!
!     Take one time step, updating q.
!     On entry, q gives
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
    use sweep_module, only: x_sweep_1st_order_gpu, x_sweep_2nd_order_gpu
    use problem_para_module, only: cc, zz, bulk, rho
#ifdef CUDA
    use cuda_module, only: threads_and_blocks
    use cuda_module, only: numBlocks, numThreads, device_id
    use memory_module, only: gpu_allocate, gpu_deallocate
    use cudafor, only: cudaMemcpy
#endif

    implicit none
    
    external rpn2, rpt2
    
    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my
    real(kind=8), intent(in) :: dx,dy,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    ! real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

! #ifdef CUDA
!     attributes(device) :: q, fm, fp, gm, gp
! #endif
    
    ! Looping scalar storage
    integer :: i,j
    real(kind=8) :: dtdx,dtdy
    real(kind=8) :: sx(mwaves, 1-mbc:mx + mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: sy(mwaves, 2-mbc:mx+mbc-1, 1-mbc:my + mbc)
    real(kind=8) :: wave_x(meqn, mwaves, 1-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: wave_y(meqn, mwaves, 2-mbc:mx+mbc-1, 1-mbc:my+mbc)
    real(kind=8), allocatable :: cflxy(:,:)
    integer :: cflmx, cflmy
#ifdef CUDA
    real(kind=8), allocatable, device :: cflxy_d(:,:)
    integer :: istat
    double precision, dimension(:,:,:), pointer, contiguous, device :: &
        q_d, fp_d, gp_d, fm_d, gm_d, sx_d
    double precision, dimension(:,:,:,:), pointer, contiguous, device :: &
        wave_x_d
    integer :: data_size
#endif

    
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

#ifdef CUDA
        call gpu_allocate( q_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)
        call gpu_allocate(fm_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)
        call gpu_allocate(fp_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)
        call gpu_allocate(sx_d, device_id, 1, mwaves, 1-mbc, mx + mbc, 2-mbc, my+mbc-1)
        call gpu_allocate(wave_x_d, device_id, 1, meqn, 1, mwaves, 1-mbc, mx+mbc, 2-mbc, my+mbc-1)
        ! call gpu_allocate(wave_x_tilde_d, device_id, 1, meqn, 1, mwaves, 1-mbc, mx+mbc, 2-mbc, my+mbc-1)
        call gpu_allocate(gm_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)
        call gpu_allocate(gp_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)

        data_size = meqn * (mx + 2*mbc) * (my + 2*mbc)

        istat = cudaMemcpy(  q_d,  q, data_size)
        istat = cudaMemcpy( fm_d, fm, data_size)
        istat = cudaMemcpy( fp_d, fp, data_size)
        istat = cudaMemcpy( sx_d, sx, mwaves*(mx+2*mbc)*(my+2*mbc-2))
        istat = cudaMemcpy( wave_x_d, wave_x, meqn*mwaves*(mx+2*mbc)*(my+2*mbc-2))
        ! istat = cudaMemcpy( wave_x_tilde_d, wave_x_d, meqn*mwaves*(mx+2*mbc)*(my+2*mbc-2))
        istat = cudaMemcpy( gm_d, gm, data_size)
        istat = cudaMemcpy( gp_d, gp, data_size)
#endif



!   -----------------------------------------------------------
!   # compute amdq and apdq
!   -----------------------------------------------------------
    call threads_and_blocks([2-mbc, 2-mbc] , [mx+mbc, my+mbc-1], numBlocks, numThreads)
    allocate(cflxy(numBlocks%x,numBlocks%y))
    allocate(cflxy_d(numBlocks%x,numBlocks%y))

    ! call x_sweep_1st_order(q, fm, fp, sx, wave_x, meqn, mwaves, mbc, mx, my, dtdx, cflgrid)

    ! call x_sweep_1st_order<<<numBlocks, numThreads, 8*numThreads%x*numThreads%y>>>&
    !     (q, fm, fp, sx, wave_x, meqn, mwaves, mbc, mx, my, dtdx, cflxy_d, cflmx, cflmy)

    call x_sweep_1st_order_gpu<<<numBlocks, numThreads, 8*numThreads%x*numThreads%y>>>&
        (q_d, fm_d, fp_d, sx_d, wave_x_d, &
         meqn, mwaves, mbc, mx, my, dtdx, cflxy_d, cflmx, cflmy, cc, zz)

    istat = cudaMemcpy(cflxy, cflxy_d, numBlocks%x*numBlocks%y)
    do i = 1,cflmx
        do j = 1,cflmy
            cflgrid = max(cflgrid, cflxy(i,j))
        enddo
    enddo

    deallocate(cflxy)
    deallocate(cflxy_d)

#ifdef CUDA
        istat = cudaMemcpy( fm,  fm_d, data_size)
        istat = cudaMemcpy( fp,  fp_d, data_size)
        istat = cudaMemcpy( gm,  gm_d, data_size)
        istat = cudaMemcpy( gp,  gp_d, data_size)
        istat = cudaMemcpy( sx,  sx_d, mwaves*(mx+2*mbc)*(my+2*mbc-2))
        istat = cudaMemcpy( wave_x, wave_x_d, meqn*mwaves*(mx+2*mbc)*(my+2*mbc-2))
#endif


    
!     -----------------------------------------------------------
!     # modify F fluxes for second order q_{xx} correction terms
!     # and solve for transverse waves
!     -----------------------------------------------------------
    call x_sweep_2nd_order(fm, fp, gm, gp, sx, wave_x, meqn, mwaves, mbc, mx, my, dtdx)

    ! call threads_and_blocks([1, 0] , [mx+1, my+1], numBlocks, numThreads)
    ! call x_sweep_2nd_order_gpu<<<numBlocks, numThreads>>>&
    !     (fm_d, fp_d, gm_d, gp_d, sx_d, wave_x_d, meqn, mwaves, mbc, mx, my, dtdx, cc, zz)

#ifdef CUDA
        call gpu_deallocate(q_d) 
        call gpu_deallocate(fp_d) 
        call gpu_deallocate(fm_d) 
        call gpu_deallocate(gp_d) 
        call gpu_deallocate(gm_d) 
        call gpu_deallocate(sx_d) 
        call gpu_deallocate(wave_x_d) 
#endif
    


!     -----------------------------------------------------------
!     # compute bmdq and bpdq
!     -----------------------------------------------------------
    call y_sweep_1st_order(q, gm, gp, sy, wave_y, meqn, mwaves, mbc, mx, my, dtdy, cflgrid)

!     -----------------------------------------------------------
!     # modify G fluxes for second order q_{yy} correction terms:
!     # and transverse waves
!     -----------------------------------------------------------
    call y_sweep_2nd_order(fm, fp, gm, gp, sy, wave_y, meqn, mwaves, mbc, mx, my, dtdy)

    ! cflgrid =  max from cfl_x and cfly

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
