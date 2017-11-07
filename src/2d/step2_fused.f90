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
    use sweep_module, only: x_sweep_1st_order, x_sweep_2nd_order, y_sweep_1st_order, y_sweep_2nd_order, x_sweep_2nd_order_simple
    use sweep_module, only: x_sweep_1st_order_gpu, x_sweep_2nd_order_gpu
    use sweep_module, only: y_sweep_1st_order_gpu, y_sweep_2nd_order_gpu
    use problem_para_module, only: cc, zz, bulk, rho
#ifdef CUDA
    use cuda_module, only: threads_and_blocks, numBlocks, numThreads, device_id
    use cuda_module, only: check_cuda_error, toString, temp_count, write_grid
    use memory_module, only: gpu_allocate, gpu_deallocate
    use cudafor, only: cudaMemcpy, cudaDeviceSynchronize
    use cudafor, only: dim3
#endif

    implicit none
    
    external rpn2, rpt2
    
    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my
    real(kind=8), intent(in) :: dx,dy,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

! #ifdef CUDA
!     attributes(device) :: q, fm, fp, gm, gp
! #endif
    
    ! Looping scalar storage
    integer :: i,j, m, mw
    real(kind=8) :: dtdx,dtdy
    real(kind=8) :: sx(mwaves, 2-mbc:mx + mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: wave_x(meqn, mwaves, 2-mbc:mx+mbc, 2-mbc:my+mbc-1)
    real(kind=8) :: sy(mwaves, 2-mbc:mx+mbc-1, 2-mbc:my + mbc)
    real(kind=8) :: wave_y(meqn, mwaves, 2-mbc:mx+mbc-1, 2-mbc:my+mbc)

    ! real(kind=8) :: cqxx(meqn,2-mbc:mx+mbc, 2-mbc:my+mbc-1)
    ! real(kind=8) :: wave_x_tilde(meqn, mwaves, 2-mbc:mx+mbc, 2-mbc:my+mbc-1)


    real(kind=8), allocatable :: cflxy(:,:)
    integer :: cflmx, cflmy
#ifdef CUDA
    real(kind=8), allocatable, device :: cflxy_d(:,:)
    integer :: istat
    double precision, dimension(:,:,:), pointer, contiguous, device :: &
        q_d, fp_d, gp_d, fm_d, gm_d, sx_d, sy_d
    double precision, dimension(:,:,:,:), pointer, contiguous, device :: &
        wave_x_d, wave_y_d
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
    ! cqxx = 0.d0
    ! wave_x_tilde = 0.d0

#ifdef CUDA
    call gpu_allocate( q_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)
    call gpu_allocate(fm_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)
    call gpu_allocate(fp_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)
    call gpu_allocate(gm_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)
    call gpu_allocate(gp_d, device_id, 1, meqn, 1-mbc,mx+mbc, 1-mbc,my+mbc)
    call gpu_allocate(sx_d, device_id, 1, mwaves, 2-mbc, mx + mbc, 2-mbc, my+mbc-1)
    call gpu_allocate(wave_x_d, device_id, 1, meqn, 1, mwaves, 2-mbc, mx+mbc, 2-mbc, my+mbc-1)
    call gpu_allocate(sy_d, device_id, 1, mwaves, 2-mbc, mx+mbc-1, 2-mbc, my+mbc)
    call gpu_allocate(wave_y_d, device_id, 1, meqn, 1, mwaves, 2-mbc, mx+mbc-1, 2-mbc, my+mbc)

    ! call gpu_allocate(cqxx_d, device_id, 1, meqn, 2-mbc,mx+mbc, 2-mbc,my+mbc-1)
    ! call gpu_allocate(wave_x_tilde_d, device_id, 1, meqn, 1, mwaves, 2-mbc, mx+mbc, 2-mbc, my+mbc-1)

    data_size = meqn * (mx + 2*mbc) * (my + 2*mbc)

    fm_d = 0.d0
    fp_d = 0.d0
    gm_d = 0.d0
    gp_d = 0.d0

    ! cqxx_d = 0.d0
    ! wave_x_tilde_d = 0.d0

    istat = cudaMemcpy(  q_d,  q, data_size)
    istat = cudaMemcpy( fm_d, fm, data_size)
    istat = cudaMemcpy( fp_d, fp, data_size)
    istat = cudaMemcpy( gm_d, gm, data_size)
    istat = cudaMemcpy( gp_d, gp, data_size)
#endif



!   -----------------------------------------------------------
!   # compute amdq and apdq
!   -----------------------------------------------------------
#ifndef gpu
    call x_sweep_1st_order(q, fm, fp, sx, wave_x, meqn, mwaves, mbc, mx, my, dtdx, cflgrid)

#else

    ! compute numBlocks and numThreads
    call threads_and_blocks([2-mbc, 0] , [mx+mbc, my+1], numBlocks, numThreads)

    cflmx = numBlocks%x
    cflmy = numBlocks%y
    allocate(cflxy(cflmx,cflmy))
    allocate(cflxy_d(cflmx,cflmy))

    call x_sweep_1st_order_gpu<<<numBlocks, numThreads, 8*numThreads%x*numThreads%y>>>&
        (q_d, fm_d, fp_d, sx_d, wave_x_d, &
         mbc, mx, my, dtdx, cflxy_d, cflmx, cflmy, cc, zz)

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

    deallocate(cflxy)
    deallocate(cflxy_d)

    istat = cudaMemcpy(fm, fm_d, data_size)
    istat = cudaMemcpy(fp, fp_d, data_size)
    istat = cudaMemcpy(gm, gm_d, data_size)
    istat = cudaMemcpy(gp, gp_d, data_size)
    istat = cudaMemcpy(sx, sx_d, mwaves*(mx + mbc-(2-mbc)+1)*( my+mbc-1-(2-mbc)+1))
    istat = cudaMemcpy(wave_x, wave_x_d, meqn*mwaves*(mx+mbc-(2-mbc)+1)*(my+mbc-1-(2-mbc)+1))

    ! output fm, fp ,sx, wave_x
    do m = 1,meqn
        open (1, file='q_after_x_1st_'//trim(toString(m,3))//'.txt', position="append")
        write(1, *) "frame: ", temp_count
        write(1, *) q(m,:,:)
        close(1)
    enddo
    do m = 1,meqn
        open (1, file='fm_after_x_1st_'//trim(toString(m,3))//'.txt', position="append")
        write(1, *) "frame: ", temp_count
        write(1, *) fm(m,:,:)
        close(1)
    enddo
    do m = 1,meqn
        open (1, file='fp_after_x_1st_'//trim(toString(m,3))//'.txt', position="append")
        write(1, *) "frame: ", temp_count
        write(1, *) fp(m,:,:)
        close(1)
    enddo
    do mw = 1,mwaves
        open (1, file='sx_after_x_1st_'//trim(toString(mw,3))//'.txt', position="append")
        write(1, *) "frame: ", temp_count
        write(1, *) sx(mw,:,:)
        close(1)
    enddo
    do mw = 1,mwaves
        do m = 1,meqn
            open (1, file='wave_x_after_x_1st_'//trim(toString(mw,3))//'_'//trim(toString(m,3))//'.txt', position="append")
            write(1, *) "frame: ", temp_count
            write(1, *) wave_x(m,mw,:,:)
            close(1)
        enddo
    enddo
    temp_count = temp_count + 1
    if (temp_count > 10) then
        stop
    endif
#endif


    
!     -----------------------------------------------------------
!     # modify F fluxes for second order q_{xx} correction terms
!     # and solve for transverse waves
!     -----------------------------------------------------------

#ifndef gpu
    ! call x_sweep_2nd_order(fm, fp, gm, gp, sx, wave_x, meqn, mwaves, mbc, mx, my, dtdx)
    call x_sweep_2nd_order_simple(fm, fp, gm, gp, sx, wave_x, meqn, mwaves, mbc, mx, my, dtdx)
    ! output fm, fp , gm, gp sx, wave_x
    ! do m = 1,meqn
    !     call write_grid( cqxx(m,:,:), 2-mbc,mx+mbc, 2-mbc, my+mbc-1, 'cqxx_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do mw = 1,mwaves
    !     do m = 1,meqn
    !         call write_grid( wave_x_tilde(m,mw,:,:), 2-mbc,mx+mbc, 2-mbc, my+mbc-1, 'wave_x_tilde_after_x_2nd_'//trim(toString(mw,3))//'_'//trim(toString(m,3))//'.txt', temp_count)
    !     enddo
    ! enddo

    ! do m = 1,meqn
    !     call write_grid( fm(m,:,:), 1-mbc,mx+mbc, 1-mbc, my+mbc, 'fm_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do m = 1,meqn
    !     call write_grid( fp(m,:,:), 1-mbc,mx+mbc, 1-mbc, my+mbc, 'fp_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do m = 1,meqn
    !     call write_grid( gm(m,:,:), 1-mbc,mx+mbc, 1-mbc, my+mbc, 'gm_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do m = 1,meqn
    !     call write_grid( gp(m,:,:), 1-mbc,mx+mbc, 1-mbc, my+mbc, 'gp_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do mw = 1,mwaves
    !     call write_grid( sx(mw,:,:), 2-mbc,mx+mbc, 2-mbc, my+mbc-1, 'sx_after_x_2nd_'//trim(toString(mw,3))//'.txt', temp_count)
    ! enddo

    ! do mw = 1,mwaves
    !     do m = 1,meqn
    !         call write_grid( wave_x(m,mw,:,:), 2-mbc,mx+mbc, 2-mbc, my+mbc-1, 'wave_x_after_x_2nd_'//trim(toString(mw,3))//'_'//trim(toString(m,3))//'.txt', temp_count)
    !     enddo
    ! enddo

    ! temp_count = temp_count + 1
    ! if (temp_count > 2) then
    !     stop
    ! endif

#else
    istat = cudaMemcpy(  q_d,  q, data_size)
    istat = cudaMemcpy( fm_d, fm, data_size)
    istat = cudaMemcpy( fp_d, fp, data_size)
    istat = cudaMemcpy( gm_d, gm, data_size)
    istat = cudaMemcpy( gp_d, gp, data_size)
    istat = cudaMemcpy( sx_d, sx, mwaves*(mx + mbc-(2-mbc)+1)*( my+mbc-1-(2-mbc)+1))
    istat = cudaMemcpy( wave_x_d, wave_x, meqn*mwaves*(mx+mbc-(2-mbc)+1)*(my+mbc-1-(2-mbc)+1))

    call threads_and_blocks([1, 0] , [mx+1, my+1], numBlocks, numThreads)
    call x_sweep_2nd_order_gpu<<<numBlocks, numThreads>>>&
        (fm_d, fp_d, gm_d, gp_d, sx_d, wave_x_d, mbc, mx, my, dtdx, cc, zz)

#ifdef DEBUG
    istat = cudaDeviceSynchronize()
    call check_cuda_error(istat)
#endif
    istat = cudaMemcpy(fm, fm_d, data_size)
    istat = cudaMemcpy(fp, fp_d, data_size)
    istat = cudaMemcpy(gm, gm_d, data_size)
    istat = cudaMemcpy(gp, gp_d, data_size)
    istat = cudaMemcpy(sx, sx_d, mwaves*(mx + mbc-(2-mbc)+1)*( my+mbc-1-(2-mbc)+1))
    istat = cudaMemcpy(wave_x, wave_x_d, meqn*mwaves*(mx+mbc-(2-mbc)+1)*(my+mbc-1-(2-mbc)+1))

    ! istat = cudaMemcpy(cqxx, cqxx_d, meqn*(mx + mbc-(2-mbc)+1)*( my+mbc-1-(2-mbc)+1))
    ! istat = cudaMemcpy(wave_x_tilde, wave_x_tilde_d, meqn*mwaves*(mx+mbc-(2-mbc)+1)*(my+mbc-1-(2-mbc)+1))

    ! output fm, fp , gm, gp sx, wave_x
    ! do m = 1,meqn
    !     call write_grid( cqxx(m,:,:), 2-mbc,mx+mbc, 2-mbc, my+mbc-1, 'cqxx_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do mw = 1,mwaves
    !     do m = 1,meqn
    !         call write_grid( wave_x_tilde(m,mw,:,:), 2-mbc,mx+mbc, 2-mbc, my+mbc-1, 'wave_x_tilde_after_x_2nd_'//trim(toString(mw,3))//'_'//trim(toString(m,3))//'.txt', temp_count)
    !     enddo
    ! enddo

    ! do m = 1,meqn
    !     call write_grid( fm(m,:,:), 1-mbc,mx+mbc, 1-mbc, my+mbc, 'fm_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do m = 1,meqn
    !     call write_grid( fp(m,:,:), 1-mbc,mx+mbc, 1-mbc, my+mbc, 'fp_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do m = 1,meqn
    !     call write_grid( gm(m,:,:), 1-mbc,mx+mbc, 1-mbc, my+mbc, 'gm_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do m = 1,meqn
    !     call write_grid( gp(m,:,:), 1-mbc,mx+mbc, 1-mbc, my+mbc, 'gp_after_x_2nd_'//trim(toString(m,3))//'.txt', temp_count)
    ! enddo

    ! do mw = 1,mwaves
    !     call write_grid( sx(mw,:,:), 2-mbc,mx+mbc, 2-mbc, my+mbc-1, 'sx_after_x_2nd_'//trim(toString(mw,3))//'.txt', temp_count)
    ! enddo

    ! do mw = 1,mwaves
    !     do m = 1,meqn
    !         call write_grid( wave_x(m,mw,:,:), 2-mbc,mx+mbc, 2-mbc, my+mbc-1, 'wave_x_after_x_2nd_'//trim(toString(mw,3))//'_'//trim(toString(m,3))//'.txt', temp_count)
    !     enddo
    ! enddo


    ! temp_count = temp_count + 1
    ! if (temp_count > 2) then
    !     stop
    ! endif

#endif




!     -----------------------------------------------------------
!     # compute bmdq and bpdq
!     -----------------------------------------------------------
#ifdef gpu
    call y_sweep_1st_order(q, gm, gp, sy, wave_y, meqn, mwaves, mbc, mx, my, dtdy, cflgrid)
#else

    istat = cudaMemcpy(  q_d,  q, data_size)
    istat = cudaMemcpy( fm_d, fm, data_size)
    istat = cudaMemcpy( fp_d, fp, data_size)
    istat = cudaMemcpy( gm_d, gm, data_size)
    istat = cudaMemcpy( gp_d, gp, data_size)

    call threads_and_blocks([0, 2-mbc] , [mx+1, my+mbc], numBlocks, numThreads)
    cflmx = numBlocks%x
    cflmy = numBlocks%y
    allocate(cflxy(cflmx,cflmy))
    allocate(cflxy_d(cflmx,cflmy))
    call y_sweep_1st_order_gpu<<<numBlocks, numThreads, 8*numThreads%x*numThreads%y>>>&
    (q_d, gm_d, gp_d, sy_d, wave_y_d, &
     mbc, mx, my, dtdy, cflxy_d, cflmx, cflmy, cc, zz)

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


    deallocate(cflxy)
    deallocate(cflxy_d)

    istat = cudaMemcpy(fm, fm_d, data_size)
    istat = cudaMemcpy(fp, fp_d, data_size)
    istat = cudaMemcpy(gm, gm_d, data_size)
    istat = cudaMemcpy(gp, gp_d, data_size)
    istat = cudaMemcpy(sy, sy_d, mwaves*(mx + mbc-1-(2-mbc)+1)*( my+mbc-(2-mbc)+1))
    istat = cudaMemcpy(wave_y, wave_y_d, meqn*mwaves*(mx+mbc-1-(2-mbc)+1)*(my+mbc-(2-mbc)+1))
#endif




!     -----------------------------------------------------------
!     # modify G fluxes for second order q_{yy} correction terms:
!     # and transverse waves
!     -----------------------------------------------------------
    call y_sweep_2nd_order(fm, fp, gm, gp, sy, wave_y, meqn, mwaves, mbc, mx, my, dtdy)

    ! call threads_and_blocks([0, 1] , [mx+1, my+1], numBlocks, numThreads)
    ! call y_sweep_2nd_order_gpu<<<numBlocks, numThreads>>>&
    !     (fm_d, fp_d, gm_d, gp_d, sy_d, wave_y_d, mbc, mx, my, dtdx, cc, zz)

#ifdef DEBUG
    istat = cudaDeviceSynchronize()
    call check_cuda_error(istat)
#endif

#ifdef CUDA
    ! istat = cudaMemcpy(fm,  fm_d, data_size)
    ! istat = cudaMemcpy(fp,  fp_d, data_size)
    ! istat = cudaMemcpy(gm,  gm_d, data_size)
    ! istat = cudaMemcpy(gp,  gp_d, data_size)
#endif



#ifdef CUDA
        call gpu_deallocate(q_d) 
        call gpu_deallocate(fp_d) 
        call gpu_deallocate(fm_d) 
        call gpu_deallocate(gp_d) 
        call gpu_deallocate(gm_d) 
        call gpu_deallocate(sx_d) 
        call gpu_deallocate(wave_x_d) 
        call gpu_deallocate(sy_d) 
        call gpu_deallocate(wave_y_d) 

        ! call gpu_deallocate(cqxx_d)
#endif
    


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
