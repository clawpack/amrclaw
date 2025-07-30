c
c --------------------------------------------------------------
c
      subroutine advanc (level,nvar,dtlevnew,vtime,naux)
c
      use amr_module
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

c     maxgr is maximum number of grids  many things are
c     dimensioned at, so this is overall. only 1d array
c     though so should suffice. problem is
c     not being able to dimension at maxthreads


c
c  ::::::::::::::; ADVANC :::::::::::::::::::::::::::::::::::::::::::
c  integrate all grids at the input  'level' by one step of its delta(t)
c  this includes:  setting the ghost cells
c                  advancing the solution on the grid
c                  adjusting fluxes for flux conservation step later
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c get start time for more detailed timing by level
      call system_clock(clock_start,clock_rate)
      call cpu_time(cpu_start)

      hx   = hxposs(level)
      delt = possk(level)
c     this is linear alg.
c     call prepgrids(listgrids,numgrids(level),level)
c

      call system_clock(clock_startBound,clock_rate)
      call cpu_time(cpu_startBound)


c     maxthreads initialized to 1 above in case no openmp
!$    maxthreads = omp_get_max_threads()

c We want to do this regardless of the threading type
!$OMP PARALLEL DO PRIVATE(j,locnew, locaux, mptr,nx,mitot,
!$OMP&                    time,levSt),
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
         mitot  = nx + 2*nghost
         locnew = node(store1,mptr)
         locaux = node(storeaux,mptr)
         time   = rnode(timemult,mptr)
c
          call bound(time,nvar,nghost,alloc(locnew),mitot,mptr,
     1               alloc(locaux),naux)

       end do
!$OMP END PARALLEL DO


      call system_clock(clock_finishBound,clock_rate)
      call cpu_time(cpu_finishBound)
      timeBound = timeBound + clock_finishBound - clock_startBound
      timeBoundCPU=timeBoundCPU+cpu_finishBound-cpu_startBound

c
c save coarse level values if there is a finer level for wave fixup

      if (level+1 .le. mxnest) then
         if (lstart(level+1) .ne. null) then
            call saveqc(level+1,nvar,naux)
         endif
      endif
c
      dtlevnew = rinfinity
      cfl_level = 0.d0    !# to keep track of max cfl seen on each level

c
      call system_clock(clock_startStepgrid,clock_rate)
      call cpu_time(cpu_startStepgrid)


!$OMP PARALLEL DO PRIVATE(j,mptr,nx,mitot)
!$OMP&            PRIVATE(mythread,dtnew)
!$OMP&            SHARED(rvol,rvoll,level,nvar,mxnest,alloc,intrat)
!$OMP&            SHARED(nghost,intratx,hx,naux,listsp)
!$OMP&            SHARED(node,rnode,dtlevnew,numgrids,listgrids)
!$OMP&            SHARED(listOfGrids,listStart,levSt)
!$OMP&            SCHEDULE (DYNAMIC,1)
!$OMP&            DEFAULT(none)
      do j = 1, numgrids(level)
         !mptr   = listgrids(j)
         levSt = listStart(level)
         mptr = listOfGrids(levSt+j-1)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         mitot  = nx + 2*nghost
c
          call par_advanc(mptr,mitot,nvar,naux,dtnew)

!$OMP CRITICAL (newdt)
          dtlevnew = dmin1(dtlevnew,dtnew)
!$OMP END CRITICAL (newdt)

      end do
!$OMP END PARALLEL DO
c

      call system_clock(clock_finish,clock_rate)
      call cpu_time(cpu_finish)
      tvoll(level) = tvoll(level) + clock_finish - clock_start
      tvollCPU(level) = tvollCPU(level) + cpu_finish - cpu_start
      timeStepgrid = timeStepgrid +clock_finish-clock_startStepgrid
      timeStepgridCPU=timeStepgridCPU+cpu_finish-cpu_startStepgrid

      cflmax = dmax1(cflmax, cfl_level)

c
      return
      end
c
c -------------------------------------------------------------
c
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
      end

c
c --------------------------------------------------------------
c
      subroutine par_advanc (mptr,mitot,nvar,naux,dtnew)
c
      use amr_module
      use gauges_module, only: update_gauges, num_gauges
      implicit double precision (a-h,o-z)


      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/

      double precision fp(nvar,mitot),fm(nvar,mitot)

c
c  :::::::::::::: PAR_ADVANC :::::::::::::::::::::::::::::::::::::::::::
c  integrate this grid. grids are done in parallel.
c  extra subr. used to allow for stack based allocation of
c  flux arrays. They are only needed temporarily. If used alloc
c  array for them it has too long a lendim, makes too big
c  a checkpoint file, and is a big critical section.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      level = node(nestlevel,mptr)
      hx    = hxposs(level)
      delt  = possk(level)
      nx    = node(ndihi,mptr) - node(ndilo,mptr) + 1
      time  = rnode(timemult,mptr)

!$    mythread = omp_get_thread_num()

      locold = node(store2, mptr)
      locnew = node(store1, mptr)

c
c  copy old soln. values into  next time step's soln. values
c  since integrator will overwrite it. only for grids not at
c  the finest level. finest level grids do not maintain copies
c  of old and new time solution values.
c
         if (level .lt. mxnest) then
             ntot   = mitot * nvar
cdir$ ivdep
             do 10 i = 1, ntot
 10            alloc(locold + i - 1) = alloc(locnew + i - 1)
         endif
c
         xlow = rnode(cornxlo,mptr) - nghost*hx

!$OMP CRITICAL(rv)
      rvol = rvol + nx
      rvoll(level) = rvoll(level) + nx
!$OMP END CRITICAL(rv)


         locaux = node(storeaux,mptr)
         call b4step1(nghost, nx, nvar, alloc(locnew),
     &                rnode(cornxlo,mptr), hx, time, dt, naux,
     &                alloc(locaux))
c
         if (node(ffluxptr,mptr) .ne. 0) then

         lenbc  = 2
            locsvf = node(ffluxptr,mptr)
            locsvq = locsvf + nvar*lenbc
            locx1d = locsvq + nvar*lenbc

         call qad(alloc(locnew),mitot,nvar,
     1             alloc(locsvf),alloc(locsvq),lenbc,
     2            intratx(level-1),hx,
     3            naux,alloc(locaux),alloc(locx1d),delt,mptr)
         endif

c        # See if the grid about to be advanced has gauge data to output.
c        # This corresponds to previous time step, but output done
c        # now to make linear interpolation easier, since grid
c        # now has boundary conditions filled in.

c     should change the way print_gauges does io - right now is critical section
c     no more,  each gauge has own array.

      if (num_gauges > 0) then
           call update_gauges(alloc(locnew:locnew+nvar*mitot),
     .                       alloc(locaux:locaux+naux*mitot),
     .                       xlow,nvar,mitot,naux,mptr)
         endif

c
         call stepgrid(alloc(locnew),fm,fp,
     2                  mitot,nghost,
     3                  delt,dtnew,hx,nvar,
     4                  xlow,time,mptr,naux,alloc(locaux))
c

      if (node(cfluxptr,mptr) .ne. 0) then
        call fluxsv(mptr,fm,fp,
     2               alloc(node(cfluxptr,mptr)),mitot,
     3               nvar,listsp(level),delt,hx)
      endif


         if (node(ffluxptr,mptr) .ne. 0) then
            lenbc = 2
            locsvf = node(ffluxptr,mptr)
         call fluxad(fm,fp,
     2               alloc(locsvf),mptr,mitot,nvar,
     4                  lenbc,intratx(level-1),
     5               nghost,delt,hx)
         endif

c
c        write(outunit,969) mythread,delt, dtnew
c969     format(" thread ",i4," updated by ",e15.7, " new dt ",e15.7)
         rnode(timemult,mptr)  = rnode(timemult,mptr)+delt
c
      return
      end
