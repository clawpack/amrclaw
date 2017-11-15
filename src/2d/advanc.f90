!> Integrate all grids at the input **level** by one step of its delta(t)
!!
!! this includes:  
!! - setting the ghost cells 
!! - advancing the solution on the grid
!! - adjusting fluxes for flux conservation step later
! --------------------------------------------------------------
!
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
    use cudafor
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
    type(grid2d) :: qs(numgrids(level)), fps(numgrids(level)), fms(numgrids(level)), &
        gps(numgrids(level)), gms(numgrids(level)) 
    type(grid2d_device) :: qs_d(numgrids(level)), fps_d(numgrids(level)), fms_d(numgrids(level)), &
        gps_d(numgrids(level)), gms_d(numgrids(level)) 
    double precision :: dtnews(numgrids(level))
    double precision, allocatable :: fm(:,:,:)
    double precision, allocatable :: fp(:,:,:)
    double precision, allocatable :: gm(:,:,:)
    double precision, allocatable :: gp(:,:,:)
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


#ifdef CUDA
    levSt = listStart(level)
    do j = 1, numgrids(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        id = j

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

        if (node(ffluxptr,mptr) .ne. 0) then
            lenbc  = 2*(nx/intratx(level-1)+ny/intraty(level-1))
            locsvf = node(ffluxptr,mptr)
            locsvq = locsvf + nvar*lenbc
            locx1d = locsvq + nvar*lenbc
            call qad(alloc(locnew),mitot,mjtot,nvar, &
                     alloc(locsvf),alloc(locsvq),lenbc, &
                     intratx(level-1),intraty(level-1),hx,hy, &
                     naux,alloc(locaux),alloc(locx1d),delt,mptr)
        endif

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

        ! allocate fluxes array and q array in SoA
        call cpu_allocate_pinned( qs(j)%dataptr, 1, mitot, 1, mjtot, 1, nvar) 
        call cpu_allocate_pinned(fms(j)%dataptr, 1, mitot, 1, mjtot, 1, nvar) 
        call cpu_allocate_pinned(fps(j)%dataptr, 1, mitot, 1, mjtot, 1, nvar) 
        call cpu_allocate_pinned(gms(j)%dataptr, 1, mitot, 1, mjtot, 1, nvar) 
        call cpu_allocate_pinned(gps(j)%dataptr, 1, mitot, 1, mjtot, 1, nvar) 

        call gpu_allocate( qs_d(j)%dataptr, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(fms_d(j)%dataptr, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(fps_d(j)%dataptr, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(gms_d(j)%dataptr, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(gps_d(j)%dataptr, device_id, 1, mitot, 1, mjtot, 1, nvar) 

        ! convert q array to SoA format
        call aos_to_soa_r2(qs(j)%dataptr, alloc(locnew), nvar, 1, mitot, 1, mjtot)

        ! copy q to GPU
        cudaResult = cudaMemcpyAsync(qs_d(j)%dataptr,  qs(j)%dataptr, nvar*mitot*mjtot, cudaMemcpyHostToDevice, get_cuda_stream(id,device_id))


        if (dimensional_split .eq. 0) then
!           # Unsplit method
            call stepgrid_soa(qs_d(j)%dataptr,fms_d(j)%dataptr,fps_d(j)%dataptr,gms_d(j)%dataptr,gps_d(j)%dataptr, &
                            mitot,mjtot,nghost, &
                            delt,dtnew,hx,hy,nvar, &
                            xlow,ylow,time,mptr,naux,alloc(locaux),numgrids(level),id)
        else if (dimensional_split .eq. 1) then
!           # Godunov splitting
            print *, "CUDA version not implemented."
            stop
            ! call stepgrid_dimSplit(alloc(locnew),fms(j)%dataptr,fps(j)%dataptr,gms(j)%dataptr,gps(j)%dataptr, &
            !              mitot,mjtot,nghost, &
            !              delt,dtnew,hx,hy,nvar, &
            !              xlow,ylow,time,mptr,naux,alloc(locaux))
        else 
!           # should never get here due to check in amr2
            write(6,*) '*** Strang splitting not supported'
            stop
        endif
        ! copy fluxes back to CPU
        cudaResult = cudaMemcpyAsync(fms(j)%dataptr,  fms_d(j)%dataptr, nvar*mitot*mjtot, cudaMemcpyDeviceToHost, get_cuda_stream(id,device_id))
        cudaResult = cudaMemcpyAsync(fps(j)%dataptr,  fps_d(j)%dataptr, nvar*mitot*mjtot, cudaMemcpyDeviceToHost, get_cuda_stream(id,device_id))
        cudaResult = cudaMemcpyAsync(gms(j)%dataptr,  gms_d(j)%dataptr, nvar*mitot*mjtot, cudaMemcpyDeviceToHost, get_cuda_stream(id,device_id))
        cudaResult = cudaMemcpyAsync(gps(j)%dataptr,  gps_d(j)%dataptr, nvar*mitot*mjtot, cudaMemcpyDeviceToHost, get_cuda_stream(id,device_id))
        cudaResult = cudaMemcpyAsync( qs(j)%dataptr,   qs_d(j)%dataptr, nvar*mitot*mjtot, cudaMemcpyDeviceToHost, get_cuda_stream(id,device_id))
    enddo

    call wait_for_all_gpu_tasks(device_id)

    do j = 1, numgrids(level)
        mptr = listOfGrids(levSt+j-1)
        nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        mitot  = nx + 2*nghost
        mjtot  = ny + 2*nghost
        locnew = node(store1, mptr)
        allocate(fm(nvar, mitot, mjtot))
        allocate(fp(nvar, mitot, mjtot))
        allocate(gm(nvar, mitot, mjtot))
        allocate(gp(nvar, mitot, mjtot))
        call soa_to_aos_r2(fm, fms(j)%dataptr, nvar, 1, mitot, 1, mjtot)
        call soa_to_aos_r2(fp, fps(j)%dataptr, nvar, 1, mitot, 1, mjtot)
        call soa_to_aos_r2(gm, gms(j)%dataptr, nvar, 1, mitot, 1, mjtot)
        call soa_to_aos_r2(gp, gps(j)%dataptr, nvar, 1, mitot, 1, mjtot)
        call soa_to_aos_r2(alloc(locnew), qs(j)%dataptr, nvar, 1, mitot, 1, mjtot)

        ! TODO: use SoA in fluxsv and fluxad as well
        if (node(cfluxptr,mptr) .ne. 0) then
            call fluxsv(mptr,fm,fp,gm,gp, &
                     alloc(node(cfluxptr,mptr)),mitot,mjtot, &
                     nvar,listsp(level),delt,hx,hy)
        endif
        if (node(ffluxptr,mptr) .ne. 0) then
            lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
            locsvf = node(ffluxptr,mptr)
            call fluxad(fm,fp,gm,gp, &
                     alloc(locsvf),mptr,mitot,mjtot,nvar, &
                        lenbc,intratx(level-1),intraty(level-1), &
                     nghost,delt,hx,hy)
        endif
        ! if (node(cfluxptr,mptr) .ne. 0) then
        !     call fluxsv(mptr,fms(j)%dataptr,fps(j)%dataptr,gms(j)%dataptr,gps(j)%dataptr, &
        !              alloc(node(cfluxptr,mptr)),mitot,mjtot, &
        !              nvar,listsp(level),delt,hx,hy)
        ! endif
        ! if (node(ffluxptr,mptr) .ne. 0) then
        !     lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
        !     locsvf = node(ffluxptr,mptr)
        !     call fluxad(fms(j)%dataptr,fps(j)%dataptr,gms(j)%dataptr,gps(j)%dataptr, &
        !              alloc(locsvf),mptr,mitot,mjtot,nvar, &
        !                 lenbc,intratx(level-1),intraty(level-1), &
        !              nghost,delt,hx,hy)
        ! endif
!
!        write(outunit,969) mythread,delt, dtnew
!969     format(" thread ",i4," updated by ",e15.7, " new dt ",e15.7)
        rnode(timemult,mptr)  = rnode(timemult,mptr)+delt

        dtlevnew = dmin1(dtlevnew,dtnew)
        deallocate(fm)
        deallocate(fp)
        deallocate(gm)
        deallocate(gp)
    enddo

    do j = 1, numgrids(level)
        call cpu_deallocated_pinned( qs(j)%dataptr) 
        call cpu_deallocated_pinned(fms(j)%dataptr) 
        call cpu_deallocated_pinned(fps(j)%dataptr) 
        call cpu_deallocated_pinned(gms(j)%dataptr) 
        call cpu_deallocated_pinned(gps(j)%dataptr) 

        call gpu_deallocate( qs_d(j)%dataptr, device_id) 
        call gpu_deallocate(fms_d(j)%dataptr, device_id) 
        call gpu_deallocate(fps_d(j)%dataptr, device_id) 
        call gpu_deallocate(gms_d(j)%dataptr, device_id) 
        call gpu_deallocate(gps_d(j)%dataptr, device_id) 
    enddo

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
        !mptr   = listgrids(j)
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
    call system_clock(clock_finish,clock_rate)
    call cpu_time(cpu_finish)
    tvoll(level) = tvoll(level) + clock_finish - clock_start
    tvollCPU(level) = tvollCPU(level) + cpu_finish - cpu_start
    timeStepgrid = timeStepgrid +clock_finish-clock_startStepgrid
    timeStepgridCPU=timeStepgridCPU+cpu_finish-cpu_startStepgrid      

    cflmax = dmax1(cflmax, cfl_level)

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


