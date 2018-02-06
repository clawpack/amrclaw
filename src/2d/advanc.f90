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
    use cuda_module, only: device_id, id_copy_cflux, toString
    use cuda_module, only: wait_for_all_gpu_tasks, wait_for_stream
    use cuda_module, only: aos_to_soa_r2, soa_to_aos_r2, get_cuda_stream
    use cuda_module, only: compute_kernel_size, numBlocks, numThreads, device_id
    use timer_module, only: take_cpu_timer, cpu_timer_start, cpu_timer_stop
    use cudafor
    use reflux_module, only: qad_cpu2, qad_gpu, &
        fluxad_fused_gpu, fluxsv_fused_gpu
    use cuda_module, only: grid_type
    use problem_para_module, only: cc, zz
#ifdef PROFILE
    use profiling_module
#endif
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
    double precision, dimension(:,:), pointer, contiguous :: cfls
    double precision, dimension(:,:), pointer, contiguous, device :: cfls_d

    type(grid_type), allocatable         :: grids(:)
    ! TODO: customized memory allocator for this?
    type(grid_type), allocatable, device :: grids_d(:)
    integer :: max_lenbc

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
    integer, parameter :: timer_fluxsv_fluxad = 12
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
#ifdef PROFILE
    call startCudaProfiler("advanc level "//toString(level),level)
#endif


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

#ifdef PROFILE
    call startCudaProfiler("bound", 24)
#endif
    ! We want to do this regardless of the threading type
    !$OMP PARALLEL DO PRIVATE(j,locnew, locaux, mptr,nx,ny,mitot, &
    !$OMP                     mjtot,time,levSt), &
    !$OMP             SHARED(level, nvar, naux, alloc, intrat, delt, &
    !$OMP                    listOfGrids,listStart,nghost, &
    !$OMP                    node,rnode,numgrids,listgrids), &
    !$OMP             SCHEDULE (dynamic,1) &
    !$OMP             DEFAULT(none)
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
#ifdef PROFILE
    call endCudaProfiler() 
#endif

    call system_clock(clock_finishBound,clock_rate)
    call cpu_time(cpu_finishBound)
    timeBound = timeBound + clock_finishBound - clock_startBound
    timeBoundCPU=timeBoundCPU+cpu_finishBound-cpu_startBound

    !
    ! save coarse level values if there is a finer level for wave fixup
#ifdef PROFILE
    call startCudaProfiler("saveqc", 24)
#endif
    if (level+1 .le. mxnest) then
        if (lstart(level+1) .ne. null) then
            call saveqc(level+1,nvar,naux)
        endif
    endif
#ifdef PROFILE
    call endCudaProfiler() 
#endif
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

    !$OMP PARALLEL DO PRIVATE(j,levSt,mptr,nx,ny,mitot,mjtot) & 
    !$OMP             PRIVATE(locold,locnew,ntot,xlow,ylow,locaux) &
    !$OMP             SHARED(numgrids,listStart,level,listOfGrids,node,ndihi,ndjhi) &
    !$OMP             SHARED(nghost,store1,store2,mxnest,alloc,rnode,cornxlo,cornylo,hx,hy) &
    !$OMP             SHARED(rvol,rvoll,storeaux,num_gauges,nvar,naux) &
    !$OMP             SHARED(timer_aos_to_soa,grid_data) &
    !$OMP             SCHEDULE (DYNAMIC,1) &
    !$OMP             DEFAULT(none)
    do j = 1, numgrids(level)
        levSt = listStart(level)
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
    
#ifdef PROFILE
        call startCudaProfiler("copy q to old", 33)
#endif
        if (level .lt. mxnest) then
            ntot   = mitot * mjtot * nvar
            do i = 1, ntot
                alloc(locold + i - 1) = alloc(locnew + i - 1)
            enddo
        endif
#ifdef PROFILE
        call endCudaProfiler()
#endif

        xlow = rnode(cornxlo,mptr) - nghost*hx
        ylow = rnode(cornylo,mptr) - nghost*hy

        !$OMP CRITICAL(rv)
        rvol = rvol + nx * ny
        rvoll(level) = rvoll(level) + nx * ny
        !$OMP END CRITICAL(rv)


        locaux = node(storeaux,mptr)
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
        call take_cpu_timer('aos_to_soa', timer_aos_to_soa)
        call cpu_timer_start(timer_aos_to_soa)
        call startCudaProfiler("aos_to_soa", 14)
#endif
        ! convert q array to SoA format
        call aos_to_soa_r2(grid_data(mptr)%ptr, alloc(locnew), nvar, 1, mitot, 1, mjtot)
#ifdef PROFILE
        call endCudaProfiler() 
        call cpu_timer_stop(timer_aos_to_soa)
#endif

    enddo
    !$OMP END PARALLEL DO

#ifdef PROFILE
    call cpu_timer_stop(timer_before_gpu_loop)

    call take_cpu_timer('gpu_loop', timer_gpu_loop)
    call cpu_timer_start(timer_gpu_loop)
#endif

#ifdef PROFILE
    call startCudaProfiler("qad and step_grid",74)
#endif
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        id = j


        xlow = rnode(cornxlo,mptr) - nghost*hx
        ylow = rnode(cornylo,mptr) - nghost*hy
        locaux = node(storeaux,mptr)

        call gpu_allocate(grid_data_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        call gpu_allocate(fms_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        call gpu_allocate(fps_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        call gpu_allocate(gms_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        call gpu_allocate(gps_d(mptr)%ptr,device_id,1,mitot,1,mjtot,1,nvar)
        call gpu_allocate(sx_d(mptr)%ptr,device_id,1,mitot-1,1,mjtot-2,1,NWAVES)
        call gpu_allocate(sy_d(mptr)%ptr,device_id,1,mitot-2,1,mjtot-1,1,NWAVES)
        call gpu_allocate(wave_x_d(mptr)%ptr,device_id,1,mitot-1,1,mjtot-2,1,NEQNS,1,NWAVES)
        call gpu_allocate(wave_y_d(mptr)%ptr,device_id,1,mitot-2,1,mjtot-1,1,NEQNS,1,NWAVES)

        ! copy q to GPU
        cudaResult = cudaMemcpyAsync(grid_data_d(mptr)%ptr, grid_data(mptr)%ptr, nvar*mitot*mjtot, cudaMemcpyHostToDevice, get_cuda_stream(id,device_id))
        
        if (associated(fflux_hh(mptr)%ptr)) then
            lenbc  = 2*(nx/intratx(level-1)+ny/intraty(level-1))
            locsvq = 1 + nvar*lenbc

            ! CPU version
            ! istat = cudaMemcpy(fflux_hh(mptr)%ptr, fflux_hd(mptr)%ptr, nvar*lenbc*2+naux*lenbc)
            ! call qad_cpu2(grid_data(mptr)%ptr,mitot,mjtot,nghost,nvar, &
            !        fflux_hh(mptr)%ptr,fflux_hh(mptr)%ptr(locsvq),lenbc, &
            !        intratx(level-1),intraty(level-1),hx,hy, &
            !        delt,mptr,cc,zz)
            ! istat = cudaMemcpy(fflux_hd(mptr)%ptr, fflux_hh(mptr)%ptr, nvar*lenbc*2+naux*lenbc)

            call compute_kernel_size(numBlocks, numThreads, &
                1,2*(nx+ny))
            call qad_gpu<<<numBlocks,numThreads,0,get_cuda_stream(id,device_id)>>>( &
                   grid_data_d(mptr)%ptr,mitot,mjtot,nghost,nvar, &
                   fflux_hd(mptr)%ptr,fflux_hd(mptr)%ptr(locsvq),lenbc, &
                   intratx(level-1),intraty(level-1),hx,hy, &
                   delt,mptr,max1d,cc,zz)
        endif


        if (dimensional_split .eq. 0) then
!           # Unsplit method
            call stepgrid_soa( &
                    grid_data_d(mptr)%ptr,fms_d(mptr)%ptr,fps_d(mptr)%ptr,gms_d(mptr)%ptr,gps_d(mptr)%ptr, &
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
    call endCudaProfiler()
#endif

#ifdef PROFILE
    call take_cpu_timer('fluxsv and fluxad', timer_fluxsv_fluxad)
    call cpu_timer_start(timer_fluxsv_fluxad)
    call startCudaProfiler("fluxsv and fluxad",12)
#endif
    allocate(grids(numgrids(level)))
    allocate(grids_d(numgrids(level)))
    max_lenbc = 0

    !$OMP PARALLEL DO PRIVATE(j,levSt,mptr,nx,ny) & 
    !$OMP             SHARED(numgrids,listStart,level,listOfGrids,node,ndihi,ndjhi,grids) &
    !$OMP             SHARED(max_lenbc,intratx,intraty) &
    !$OMP             SHARED(fms_d,fps_d,gms_d,gps_d) &
    !$OMP             SCHEDULE (DYNAMIC,1) &
    !$OMP             DEFAULT(none)
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr = listOfGrids(levSt+j-1)
        nx   = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny   = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        grids(j)%fm => fms_d(mptr)%ptr
        grids(j)%fp => fps_d(mptr)%ptr
        grids(j)%gm => gms_d(mptr)%ptr
        grids(j)%gp => gps_d(mptr)%ptr
        grids(j)%mptr = mptr
        grids(j)%nx   = nx 
        grids(j)%ny   = ny 
        if (level > 1) then
            !$OMP CRITICAL (max_lenbc)
            max_lenbc = max(max_lenbc, 2*(nx/intratx(level-1)+ny/intraty(level-1)))
            !$OMP END CRITICAL (max_lenbc)
        endif
    enddo
    !$OMP END PARALLEL DO

    grids_d = grids
    ! one kernel launch to do fluxsv for all grids at this level
    ! we don't do this for then fineset level
    if (level < lfine) then
        call compute_kernel_size(numBlocks, numThreads, &
            1,listsp(level),1,numgrids(level))
        call fluxsv_fused_gpu<<<numBlocks,numThreads>>>( &
                 grids_d, cflux_dd, fflux_dd, &
                 nghost, numgrids(level), nvar,listsp(level),delt,hx,hy)
    endif
    
    ! one kernel launch to do fluxad for all grids at this level
    ! we don't do this for the coarsest level
    if (level > 1) then
        call compute_kernel_size(numBlocks, numThreads, &
            1,max_lenbc,1,numgrids(level))
        call fluxad_fused_gpu<<<numBlocks,numThreads>>>( &
                grids_d, fflux_dd,&
                nghost, numgrids(level), intratx(level-1), intraty(level-1), &
                delt, hx, hy)
        call wait_for_all_gpu_tasks(device_id)
        do j = 1, numgrids(level)
            levSt = listStart(level)
            mptr = listOfGrids(levSt+j-1)
            nx   = node(ndihi,mptr) - node(ndilo,mptr) + 1
            ny   = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
            lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
            istat = cudaMemcpy(fflux_hh(mptr)%ptr, fflux_hd(mptr)%ptr, nvar*lenbc*2+naux*lenbc)
        enddo
    endif
    call wait_for_all_gpu_tasks(device_id)
    deallocate(grids)
    deallocate(grids_d)
#ifdef PROFILE
    call endCudaProfiler()
    call cpu_timer_stop(timer_fluxsv_fluxad)
#endif

#ifdef PROFILE
    call cpu_timer_stop(timer_gpu_loop)

    call startCudaProfiler('soa_to_aos',99)
    call take_cpu_timer('soa_to_aos', timer_soa_to_aos)
    call cpu_timer_start(timer_soa_to_aos)
#endif

    !$OMP PARALLEL DO PRIVATE(j,levSt,mptr,nx,ny,mitot,mjtot) & 
    !$OMP             SHARED(numgrids,listStart,level,listOfGrids,node,ndihi,ndjhi) &
    !$OMP             SHARED(nghost,locnew,store1,alloc) &
    !$OMP             SHARED(grid_data,nvar,rnode,timebult,delt,device_id) &
    !$OMP             SHARED(grid_data_d,fms_d,fps_d,gms_d,gps_d,sx_d,sy_d,wave_x_d,wave_y_d) &
    !$OMP             SCHEDULE (DYNAMIC,1) &
    !$OMP             DEFAULT(none)
    do j = 1, numgrids(level)
        levSt = listStart(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        locnew = node(store1, mptr)

        ! TODO: use callback function to let this get executed right after grid mptr
        ! is done in its cuda stream
        call soa_to_aos_r2(alloc(locnew), grid_data(mptr)%ptr, nvar, 1, mitot, 1, mjtot)

        rnode(timemult,mptr)  = rnode(timemult,mptr)+delt

        call gpu_deallocate(grid_data_d(mptr)%ptr,device_id)
        call gpu_deallocate(fms_d(mptr)%ptr,device_id)
        call gpu_deallocate(fps_d(mptr)%ptr,device_id)
        call gpu_deallocate(gms_d(mptr)%ptr,device_id)
        call gpu_deallocate(gps_d(mptr)%ptr,device_id)
        call gpu_deallocate(sx_d(mptr)%ptr,device_id)
        call gpu_deallocate(sy_d(mptr)%ptr,device_id)
        call gpu_deallocate(wave_x_d(mptr)%ptr,device_id)
        call gpu_deallocate(wave_y_d(mptr)%ptr,device_id)
    enddo
    !$OMP END PARALLEL DO
#ifdef PROFILE
    call cpu_timer_stop(timer_soa_to_aos)
    call endCudaProfiler()
#endif




#else
    !$OMP PARALLEL DO PRIVATE(j,mptr,nx,ny,mitot,mjtot) & 
    !$OMP             PRIVATE(mythread,dtnew) &
    !$OMP             SHARED(rvol,rvoll,level,nvar,mxnest,alloc,intrat) &
    !$OMP             SHARED(nghost,intratx,intraty,hx,hy,naux,listsp) &
    !$OMP             SHARED(node,rnode,dtlevnew,numgrids,listgrids) &
    !$OMP             SHARED(listOfGrids,listStart,levSt) &
    !$OMP             SCHEDULE (DYNAMIC,1) &
    !$OMP             DEFAULT(none)
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
    call startCudaProfiler('CFL reduction', 34)
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
    call endCudaProfiler()
#endif

#else
    cflmax = dmax1(cflmax, cfl_level)
#endif

    !
#ifdef PROFILE
    call endCudaProfiler()
#endif
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

