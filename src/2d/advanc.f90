!> Integrate all grids at the input **level** by one step of its delta(t)
!!
!! this includes:  
!! - setting the ghost cells 
!! - advancing the solution on the grid
!! - adjusting fluxes for flux conservation step later
! --------------------------------------------------------------
!
#include "amr_macros.H"

subroutine advanc(level,nvar,dtlevnew,vtime,naux)
    use amr_module 
    use parallel_advanc_module
#ifdef CUDA
    use gauges_module, only: update_gauges, num_gauges
    use memory_module, only: cpu_allocate_pinned, cpu_deallocated_pinned, &
        gpu_allocate, gpu_deallocate
    use cuda_module, only: grid2d, grid2d_device, device_id
    use cuda_module, only: wait_for_all_gpu_tasks
    use cuda_module, only: aos_to_soa_r2, soa_to_aos_r2, get_cuda_stream
    use cuda_module, only: compute_kernel_size, numBlocks, numThreads, device_id
    use timer_module, only: take_cpu_timer, cpu_timer_start, cpu_timer_stop
    use cudafor
    use sweep_module, only: fluxad_gpu, fluxsv_gpu
#endif
    implicit double precision (a-h,o-z)


    logical    vtime
    integer omp_get_thread_num, omp_get_max_threads
    integer mythread/0/, maxthreads/1/
    integer listgrids(numgrids(level))
    integer clock_start, clock_finish, clock_rate
    integer clock_startStepgrid,clock_startBound,clock_finishBound
    real(kind=8) cpu_start, cpu_finish
    real(kind=8) cpu_startBound, cpu_finishBound
    real(kind=8) cpu_startStepgrid, cpu_finishStepgrid

#ifdef CUDA
    integer :: locold, locnew, locaux
    integer :: i,j, id
    integer :: cudaResult
    double precision :: xlow, ylow
    double precision :: cfl_local
    type(grid2d_device) :: fps_d(numgrids(level)), fms_d(numgrids(level)), &
        gps_d(numgrids(level)), gms_d(numgrids(level)) 
    double precision, dimension(:,:), pointer, contiguous :: cfls
    double precision, dimension(:,:), pointer, contiguous, device :: cfls_d

#ifdef PROFILE
    integer, parameter :: timer_stepgrid = 1
    integer, parameter :: timer_gpu_loop = 2
    integer, parameter :: timer_before_gpu_loop = 3
    integer, parameter :: timer_post = 4
    integer, parameter :: timer_cfl = 5
    integer, parameter :: timer_aos_to_soa = 6
    integer, parameter :: timer_soa_to_aos = 7
    integer, parameter :: timer_init_cfls = 8
    integer, parameter :: timer_qad = 9
    integer, parameter :: timer_allocate = 12
#endif
#endif

    !     maxgr is maximum number of grids  many things are
    !     dimensioned at, so this is overall. only 1d array
    !     though so should suffice. problem is
    !     not being able to dimension at maxthreads


    !
    !  ::::::::::::::: ADVANC :::::::::::::::::::::::::::::::::::::::::::
    !  integrate all grids at the input  'level' by one step of its delta(t)
    !  this includes:  setting the ghost cells 
    !                  advancing the solution on the grid
    !                  adjusting fluxes for flux conservation step later
    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    ! get start time for more detailed timing by level
    call system_clock(clock_start,clock_rate)
    call cpu_time(cpu_start)


    hx   = hxposs(level)
    hy   = hyposs(level)
    delt = possk(level)
    !     this is linear alg.
    !     call prepgrids(listgrids,numgrids(level),level)
    !

    call system_clock(clock_startBound,clock_rate)
    call cpu_time(cpu_startBound)


    !     maxthreads initialized to 1 above in case no openmp
    !$    maxthreads = omp_get_max_threads()

    ! We want to do this regardless of the threading type
    !$OMP PARALLEL DO PRIVATE(j,locnew, locaux, mptr,nx,ny,mitot,
    !$OMP&                    mjtot,time,levSt),
    !$OMP&            SHARED(level, nvar, naux, alloc, intrat, delt,
    !$OMP&                   listOfGrids,listStart,nghost,
    !$OMP&                   node,rnode,numgrids,listgrids),
    !$OMP&            SCHEDULE (dynamic,1)
    !$OMP&            DEFAULT(none)
    do j = 1, numgrids(level)
        !mptr   = listgrids(j)
        levSt = listStart(level)
        mptr   = listOfGrids(levSt+j-1)
        !write(*,*)"old ",listgrids(j)," new",mptr
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        locnew = node(store1,mptr)
        locaux = node(storeaux,mptr)
        time   = rnode(timemult,mptr)
        !     
        call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr, alloc(locaux),naux)

    end do
    !$OMP END PARALLEL DO
    call system_clock(clock_finishBound,clock_rate)
    call cpu_time(cpu_finishBound)
    timeBound = timeBound + clock_finishBound - clock_startBound
    timeBoundCPU=timeBoundCPU+cpu_finishBound-cpu_startBound

    !
    ! save coarse level values if there is a finer level for wave fixup
    if (level+1 .le. mxnest) then
        if (lstart(level+1) .ne. null) then
            call saveqc(level+1,nvar,naux)
        endif
    endif
    !
    dtlevnew = rinfinity
    cfl_level = 0.d0    !# to keep track of max cfl seen on each level

    ! 
    call system_clock(clock_startStepgrid,clock_rate)
    call cpu_time(cpu_startStepgrid)

#ifdef PROFILE
    call take_cpu_timer('stepgrid', timer_stepgrid)
    call cpu_timer_start(timer_stepgrid)
#endif



#ifdef CUDA

#ifdef PROFILE
    call take_cpu_timer('Initialize cfls', timer_init_cfls)
    call cpu_timer_start(timer_init_cfls)
#endif

    call cpu_allocate_pinned(cfls,1,numgrids(level),1,2)
    call gpu_allocate(cfls_d,device_id,1,numgrids(level),1,2)

    ! TODO: merge this with something else
    cfls_d = 0.d0

#ifdef PROFILE
    call cpu_timer_stop(timer_init_cfls)
#endif


#ifdef PROFILE
    call take_cpu_timer('pre-process before gpu loop', timer_before_gpu_loop)
    call cpu_timer_start(timer_before_gpu_loop)
#endif
    levSt = listStart(level)
    do j = 1, numgrids(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost

        !  copy old soln. values into  next time step's soln. values
        !  since integrator will overwrite it. only for grids not at
        !  the finest level. finest level grids do not maintain copies
        !  of old and new time solution values.

        locold = node(store2, mptr)
        locnew = node(store1, mptr)
    
        if (level .lt. mxnest) then
            ntot   = mitot * mjtot * nvar
            do i = 1, ntot
                alloc(locold + i - 1) = alloc(locnew + i - 1)
            enddo
        endif
        xlow = rnode(cornxlo,mptr) - nghost*hx
        ylow = rnode(cornylo,mptr) - nghost*hy

        rvol = rvol + nx * ny
        rvoll(level) = rvoll(level) + nx * ny


        locaux = node(storeaux,mptr)

#ifdef PROFILE
        call take_cpu_timer('qad', timer_qad)
        call cpu_timer_start(timer_qad)
#endif
        if (associated(fflux(mptr)%ptr)) then
            lenbc  = 2*(nx/intratx(level-1)+ny/intraty(level-1))
            locsvq = 1 + nvar*lenbc
            locx1d = locsvq + nvar*lenbc
            call qad(alloc(locnew),mitot,mjtot,nvar, &
                     fflux(mptr)%ptr,fflux(mptr)%ptr(locsvq),lenbc, &
                     intratx(level-1),intraty(level-1),hx,hy, &
                     naux,alloc(locaux),fflux(mptr)%ptr(locx1d),delt,mptr)
        endif
#ifdef PROFILE
        call cpu_timer_stop(timer_qad)
#endif

        !        # See if the grid about to be advanced has gauge data to output.
        !        # This corresponds to previous time step, but output done
        !        # now to make linear interpolation easier, since grid
        !        # now has boundary conditions filled in.

        !     should change the way print_gauges does io - right now is critical section
        !     no more,  each gauge has own array.

        if (num_gauges > 0) then
            call update_gauges(alloc(locnew:locnew+nvar*mitot*mjtot), &
                               alloc(locaux:locaux+nvar*mitot*mjtot), &
                               xlow,ylow,nvar,mitot,mjtot,naux,mptr)
        endif

#ifdef PROFILE
        call take_cpu_timer('allocate q and fluxes', timer_allocate)
        call cpu_timer_start(timer_allocate)
#endif
        ! allocate fluxes array and q array in SoA
        call gpu_allocate(fms_d(j)%dataptr, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(fps_d(j)%dataptr, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(gms_d(j)%dataptr, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(gps_d(j)%dataptr, device_id, 1, mitot, 1, mjtot, 1, nvar) 
#ifdef PROFILE
        call cpu_timer_stop(timer_allocate)
#endif

#ifdef PROFILE
        call take_cpu_timer('aos_to_soa', timer_aos_to_soa)
        call cpu_timer_start(timer_aos_to_soa)
#endif

        ! convert q array to SoA format
        call aos_to_soa_r2(grid_data(mptr)%ptr, alloc(locnew), nvar, 1, mitot, 1, mjtot)

#ifdef PROFILE
        call cpu_timer_stop(timer_aos_to_soa)
#endif
    enddo

#ifdef PROFILE
    call cpu_timer_stop(timer_before_gpu_loop)

    call take_cpu_timer('gpu_loop', timer_gpu_loop)
    call cpu_timer_start(timer_gpu_loop)
#endif

    do j = 1, numgrids(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        id = j

        locold = node(store2, mptr)
        locnew = node(store1, mptr)

        xlow = rnode(cornxlo,mptr) - nghost*hx
        ylow = rnode(cornylo,mptr) - nghost*hy
        locaux = node(storeaux,mptr)

        ! copy q to GPU
        cudaResult = cudaMemcpyAsync(grid_data_d(mptr)%ptr, grid_data(mptr)%ptr, nvar*mitot*mjtot, cudaMemcpyHostToDevice, get_cuda_stream(id,device_id))


        if (dimensional_split .eq. 0) then
!           # Unsplit method
            call stepgrid_soa(grid_data_d(mptr)%ptr,fms_d(j)%dataptr,fps_d(j)%dataptr,gms_d(j)%dataptr,gps_d(j)%dataptr, &
                            mitot,mjtot,nghost, &
                            delt,hx,hy,nvar, &
                            xlow,ylow,time,mptr,naux,alloc(locaux),& 
                            numgrids(level),id,cfls_d)
        else if (dimensional_split .eq. 1) then
!           # Godunov splitting
            print *, "CUDA version not implemented."
            stop
        else 
!           # should never get here due to check in amr2
            write(6,*) '*** Strang splitting not supported'
            stop
        endif
        cudaResult = cudaMemcpyAsync(grid_data(mptr)%ptr, grid_data_d(mptr)%ptr, nvar*mitot*mjtot, cudaMemcpyDeviceToHost, get_cuda_stream(id,device_id))
    enddo

    call wait_for_all_gpu_tasks(device_id)

#ifdef PROFILE
    call cpu_timer_stop(timer_gpu_loop)

    call take_cpu_timer('post-process after gpu loop', timer_post)
    call cpu_timer_start(timer_post)
#endif
    do j = 1, numgrids(level)
        id = j
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        locnew = node(store1, mptr)

#ifdef PROFILE
        call take_cpu_timer('soa_to_aos', timer_soa_to_aos)
        call cpu_timer_start(timer_soa_to_aos)
#endif

        call soa_to_aos_r2(alloc(locnew), grid_data(mptr)%ptr, nvar, 1, mitot, 1, mjtot)

#ifdef PROFILE
        call cpu_timer_stop(timer_soa_to_aos)
#endif

        if (associated(cflux(mptr)%ptr)) then

            ! TODO: this line does some redundant work
            cudaResult = cudaMemcpy(cflux_d(mptr)%ptr, cflux(mptr)%ptr, 5*listsp(level))
            call compute_kernel_size(numBlocks,numThreads,1,listsp(level))

            call fluxsv_gpu<<<numBlocks,numThreads>>>(mptr, &
                     fms_d(j)%dataptr,fps_d(j)%dataptr,gms_d(j)%dataptr,gps_d(j)%dataptr, &
                     cflux_d(mptr)%ptr, &
                     fflux, &
                     mitot,mjtot,nvar,listsp(level),delt,hx,hy)
        endif
    
        if (associated(fflux(mptr)%ptr)) then
            lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
            call compute_kernel_size(numBlocks, numThreads,1,lenbc)
            call fluxad_gpu<<<numBlocks,numThreads>>>(fms_d(j)%dataptr, fps_d(j)%dataptr, gms_d(j)%dataptr, gps_d(j)%dataptr, &
                nghost, nx, ny, lenbc, &
                intratx(level-1), intraty(level-1), &
                fflux(mptr)%ptr, delt, hx, hy)
        endif
        rnode(timemult,mptr)  = rnode(timemult,mptr)+delt


    enddo
    ! TODO: remove this barrier
    call wait_for_all_gpu_tasks(device_id)

    do j = 1, numgrids(level)
        call gpu_deallocate(fms_d(j)%dataptr, device_id) 
        call gpu_deallocate(fps_d(j)%dataptr, device_id) 
        call gpu_deallocate(gms_d(j)%dataptr, device_id) 
        call gpu_deallocate(gps_d(j)%dataptr, device_id) 
    enddo

#ifdef PROFILE
    call cpu_timer_stop(timer_post)
#endif

#else
    !$OMP PARALLEL DO PRIVATE(j,mptr,nx,ny,mitot,mjtot)  
    !$OMP&            PRIVATE(mythread,dtnew)
    !$OMP&            SHARED(rvol,rvoll,level,nvar,mxnest,alloc,intrat)
    !$OMP&            SHARED(nghost,intratx,intraty,hx,hy,naux,listsp)
    !$OMP&            SHARED(node,rnode,dtlevnew,numgrids,listgrids)
    !$OMP&            SHARED(listOfGrids,listStart,levSt)
    !$OMP&            SCHEDULE (DYNAMIC,1)
    !$OMP&            DEFAULT(none)
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        !
        call par_advanc(mptr,mitot,mjtot,nvar,naux,dtnew)
        !$OMP CRITICAL (newdt)
        dtlevnew = dmin1(dtlevnew,dtnew)
        !$OMP END CRITICAL (newdt)    
    end do
    !$OMP END PARALLEL DO

#endif
    !
#ifdef PROFILE
    call cpu_timer_stop(timer_stepgrid)
#endif
    call system_clock(clock_finish,clock_rate)
    call cpu_time(cpu_finish)
    tvoll(level) = tvoll(level) + clock_finish - clock_start
    tvollCPU(level) = tvollCPU(level) + cpu_finish - cpu_start
    timeStepgrid = timeStepgrid +clock_finish-clock_startStepgrid
    timeStepgridCPU=timeStepgridCPU+cpu_finish-cpu_startStepgrid      

#ifdef CUDA

#ifdef PROFILE
    call take_cpu_timer('CFL reduction', timer_cfl)
    call cpu_timer_start(timer_cfl)
#endif

    ! reduction to get cflmax and dtlevnew
    cudaResult = cudaMemcpy(cfls, cfls_d, numgrids(level)*2)
    do j = 1,numgrids(level)
        cfl_local = max(cfls(j,1),cfls(j,2))
        if (cfl_local .gt. cflv1) then
            write(*,810) cfl_local
            write(outunit,810) cfl_local, cflv1
      810   format('*** WARNING *** Courant number  =', d12.4, &
          '  is larger than input cfl_max = ', d12.4)
        endif
        cfl_level = dmax1(cfl_level, cfl_local)
    enddo
    dtlevnew = delt*cfl/cfl_level
    cflmax = dmax1(cflmax, cfl_level)
    call cpu_deallocated_pinned(cfls)
    call gpu_deallocate(cfls_d,device_id)

#ifdef PROFILE
    call cpu_timer_stop(timer_cfl)
#endif

#else
    cflmax = dmax1(cflmax, cfl_level)
#endif

    !
    return
end subroutine advanc
    !
    ! -------------------------------------------------------------
    !
subroutine prepgrids(listgrids, num, level)

    use amr_module
    implicit double precision (a-h,o-z)
    integer listgrids(num)

    mptr = lstart(level)
    do j = 1, num
        listgrids(j) = mptr
        mptr = node(levelptr, mptr)
    end do

    if (mptr .ne. 0) then
        write(*,*)" Error in routine setting up grid array "
        stop
    endif

    return
end subroutine prepgrids

