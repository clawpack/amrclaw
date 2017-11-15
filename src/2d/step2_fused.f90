!> Compute all fluxes at cell edges 
!! \param q[in] solution array for computing fluxes. It is not changed in this subroutine
!! \param fm[out] fluxes on the left side of each vertical edge
!! \param fp[out] fluxes on the right side of each vertical edge
!! \param gm[out] fluxes on the lower side of each horizontal edge
!! \param gp[out] fluxes on the upper side of each horizontal edge
subroutine step2_fused(maxm,meqn,maux,mbc,mx,my,q,dx,dy,dt,cflgrid,fm,fp,gm,gp,rpn2,rpt2,id)
!
!     clawpack routine ...  modified for AMRCLAW
!
!     Take one time step, updating q.
!     On entry, q gives
!        initial data for this step
!        and is unchanged
!    
!     fm, fp are fluxes to left and right of single cell edge
!     gm, gp are fluxes to lower and upper of single cell edge

!     q, fm, fp, gm, gp are all in Structure of Array (SoA) format
    
    use amr_module
    use parallel_advanc_module, only: dtcom, dxcom, dycom, icom, jcom
    use sweep_module, only: x_sweep_1st_order, x_sweep_2nd_order, y_sweep_1st_order, y_sweep_2nd_order, x_sweep_2nd_order_simple
    use sweep_module, only: x_sweep_1st_order_gpu, x_sweep_2nd_order_gpu
    use sweep_module, only: y_sweep_1st_order_gpu, y_sweep_2nd_order_gpu
    use problem_para_module, only: cc, zz, bulk, rho
#ifdef CUDA
    use cuda_module, only: threads_and_blocks, numBlocks, numThreads, device_id
    use cuda_module, only: check_cuda_error, write_grid
    use cuda_module, only: get_cuda_stream
    use memory_module, only: gpu_allocate, gpu_deallocate
    use memory_module, only: cpu_allocate_pinned, cpu_deallocated_pinned
    use cudafor, only: cudaMemcpyAsync, cudaDeviceSynchronize, cudaMemcpy
    use cudafor, only: dim3
#endif

    implicit none
    
    external rpn2, rpt2
    
    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my
    integer, intent(in) :: id
    real(kind=8), intent(in) :: dx,dy,dt
    real(kind=8), intent(inout) :: cflgrid

    real(kind=8), intent(in)    ::  q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
    real(kind=8), intent(inout) :: fm(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
    real(kind=8), intent(inout) :: fp(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
    real(kind=8), intent(inout) :: gm(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)
    real(kind=8), intent(inout) :: gp(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)

#ifdef CUDA
    attributes(device) :: q, fm, fp, gm, gp
#endif
    
    ! Looping scalar storage
    integer :: i,j, m, mw
    real(kind=8) :: dtdx,dtdy

#ifdef CUDA
    double precision, dimension(:,:), pointer, contiguous :: &
        cflxy
    double precision, dimension(:,:), pointer, contiguous, device :: &
        cflxy_d
    double precision, dimension(:,:,:), pointer, contiguous, device :: &
        sx_d, sy_d
    double precision, dimension(:,:,:,:), pointer, contiguous, device :: &
        wave_x_d, wave_y_d
    integer :: data_size
    integer :: istat
    integer :: cflmx, cflmy
#endif

    
    ! Store mesh parameters in common block
    dxcom = dx
    dycom = dy
    dtcom = dt
    
    cflgrid = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy

#ifdef CUDA
    call gpu_allocate(    sx_d, device_id, 2-mbc, mx+mbc, 2-mbc, my+mbc-1, 1, mwaves)
    call gpu_allocate(wave_x_d, device_id, 2-mbc, mx+mbc, 2-mbc, my+mbc-1, 1, meqn, 1, mwaves)
    call gpu_allocate(    sy_d, device_id, 2-mbc, mx+mbc-1, 2-mbc, my+mbc, 1, mwaves)
    call gpu_allocate(wave_y_d, device_id, 2-mbc, mx+mbc-1, 2-mbc, my+mbc, 1, meqn, 1, mwaves)

    data_size = meqn * (mx + 2*mbc) * (my + 2*mbc)

    ! Initialize fluxes array to zero
    call init_fluxes(fm, fp, gm, gp, meqn, 1-mbc, mx+mbc, 1-mbc, my+mbc, id)
#else
    fm = 0.d0
    fp = 0.d0
    gm = 0.d0
    gp = 0.d0

#endif



!   -----------------------------------------------------------
!   # compute amdq and apdq
!   -----------------------------------------------------------
#ifdef gpu
    call x_sweep_1st_order(q, fm, fp, sx, wave_x, meqn, mwaves, mbc, mx, my, dtdx, cflgrid)

#else

    ! compute numBlocks and numThreads
    call threads_and_blocks([2-mbc, 0] , [mx+mbc, my+1], numBlocks, numThreads)

    cflmx = numBlocks%x
    cflmy = numBlocks%y
    call cpu_allocate_pinned(cflxy, 1, cflmx,1, cflmy)
    call gpu_allocate(cflxy_d, device_id, 1, cflmx,1, cflmy)

    call x_sweep_1st_order_gpu<<<numBlocks, numThreads, 8*numThreads%x*numThreads%y, get_cuda_stream(id, device_id)>>>&
        (q, fm, fp, sx_d, wave_x_d, &
         mbc, mx, my, dtdx, cflxy_d, cc, zz)

#ifdef DEBUG
    istat = cudaDeviceSynchronize()
    call check_cuda_error(istat)
#endif


    ! TODO: replace this reduction with library version
    ! TODO: move this to advanc
    istat = cudaMemcpy(cflxy, cflxy_d, cflmx*cflmy)

    do i = 1,cflmx
        do j = 1,cflmy
            cflgrid = max(cflgrid, cflxy(i,j))
        enddo
    enddo

    call cpu_deallocated_pinned(cflxy)
    call gpu_deallocate(cflxy_d, device_id)

#endif


    
!     -----------------------------------------------------------
!     # modify F fluxes for second order q_{xx} correction terms
!     # and solve for transverse waves
!     -----------------------------------------------------------

#ifdef gpu
    ! call x_sweep_2nd_order(fm, fp, gm, gp, sx, wave_x, meqn, mwaves, mbc, mx, my, dtdx)
    call x_sweep_2nd_order_simple(fm, fp, gm, gp, sx, wave_x, meqn, mwaves, mbc, mx, my, dtdx)

#else
    call threads_and_blocks([1, 0] , [mx+1, my+1], numBlocks, numThreads)
    call x_sweep_2nd_order_gpu<<<numBlocks, numThreads,0,get_cuda_stream(id,device_id)>>>&
        (fm, fp, gm, gp, sx_d, wave_x_d, mbc, mx, my, dtdx, cc, zz)

#ifdef DEBUG
    istat = cudaDeviceSynchronize()
    call check_cuda_error(istat)
#endif
#endif




!     -----------------------------------------------------------
!     # compute bmdq and bpdq
!     -----------------------------------------------------------
#ifdef gpu
    call y_sweep_1st_order(q, gm, gp, sy, wave_y, meqn, mwaves, mbc, mx, my, dtdy, cflgrid)
#else

    call threads_and_blocks([0, 2-mbc] , [mx+1, my+mbc], numBlocks, numThreads)
    cflmx = numBlocks%x
    cflmy = numBlocks%y
    call cpu_allocate_pinned(cflxy, 1, cflmx,1, cflmy)
    call gpu_allocate(cflxy_d, device_id, 1, cflmx,1, cflmy)
    call y_sweep_1st_order_gpu<<<numBlocks, numThreads, 8*numThreads%x*numThreads%y,get_cuda_stream(id,device_id)>>>&
    (q, gm, gp, sy_d, wave_y_d, &
     mbc, mx, my, dtdy, cflxy_d, cc, zz)

#ifdef DEBUG
    istat = cudaDeviceSynchronize()
    call check_cuda_error(istat)
#endif

    istat = cudaMemcpy(cflxy, cflxy_d, cflmx*cflmy)

    do i = 1,cflmx
        do j = 1,cflmy
            cflgrid = max(cflgrid, cflxy(i,j))
        enddo
    enddo

    call cpu_deallocated_pinned(cflxy)
    call gpu_deallocate(cflxy_d, device_id)

#endif




!     -----------------------------------------------------------
!     # modify G fluxes for second order q_{yy} correction terms:
!     # and transverse waves
!     -----------------------------------------------------------
#ifdef gpu
    call y_sweep_2nd_order(fm, fp, gm, gp, sy, wave_y, meqn, mwaves, mbc, mx, my, dtdy)

#else

    call threads_and_blocks([0, 1] , [mx+1, my+1], numBlocks, numThreads)
    call y_sweep_2nd_order_gpu<<<numBlocks, numThreads,0,get_cuda_stream(id,device_id)>>>&
        (fm, fp, gm, gp, sy_d, wave_y_d, mbc, mx, my, dtdx, cc, zz)

#ifdef DEBUG
    istat = cudaDeviceSynchronize()
    call check_cuda_error(istat)
#endif

#endif



#ifdef CUDA
        call gpu_deallocate(    sx_d, device_id) 
        call gpu_deallocate(wave_x_d, device_id) 
        call gpu_deallocate(    sy_d, device_id) 
        call gpu_deallocate(wave_y_d, device_id) 
#endif
end subroutine step2_fused

! Set fluxes array to zero
! assume input fluxes are in SoA format
subroutine init_fluxes(fm, fp, gm, gp, nvar, xlo, xhi, ylo, yhi &
#ifdef CUDA
        , id &
#endif
        )

#ifdef CUDA
    use cuda_module, only: device_id
    use cuda_module, only: get_cuda_stream
#endif
    implicit none
    integer, intent(in) :: nvar, xlo, xhi, ylo, yhi
#ifdef CUDA
    ! id of this grid on current level
    integer, intent(in) :: id
#endif
    double precision, intent(inout) :: fm(xlo:xhi, ylo:yhi, 1:nvar)
    double precision, intent(inout) :: fp(xlo:xhi, ylo:yhi, 1:nvar)
    double precision, intent(inout) :: gm(xlo:xhi, ylo:yhi, 1:nvar)
    double precision, intent(inout) :: gp(xlo:xhi, ylo:yhi, 1:nvar)
    integer :: i,j,m
#ifdef CUDA
    attributes(device) :: fm, fp, gm, gp
#endif

    !$cuf kernel do(3) <<<*,*,0, get_cuda_stream(id,device_id)>>> 
    do m = 1,nvar
        do j = ylo,yhi
            do i = xlo,xhi
                fm(i,j,m) = 0.d0
                fp(i,j,m) = 0.d0
                gm(i,j,m) = 0.d0
                gp(i,j,m) = 0.d0
            enddo
        enddo
    enddo


end subroutine init_fluxes

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
