!> Integrate all grids at the input **level** by one step of its delta(t)
!!
!! this includes:  
!! - setting the ghost cells 
!! - advancing the solution on the grid
!! - adjusting fluxes for flux conservation step later
! --------------------------------------------------------------
!
      subroutine advanc (level,nvar,dtlevnew,vtime,naux)
!
      use amr_module 
      use parallel_advanc_module
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

!     maxgr is maximum number of grids  many things are
!     dimensioned at, so this is overall. only 1d array
!     though so should suffice. problem is
!     not being able to dimension at maxthreads


!
!  ::::::::::::::; ADVANC :::::::::::::::::::::::::::::::::::::::::::
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
          call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr,
     1               alloc(locaux),naux)

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
      end
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
      end
