
c --------------------------------------------------------------
c
      subroutine advanc (level,nvar,dtlevnew,vtime,naux)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      logical    vtime
      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/
      integer listgrids(numgrids(level))

c     maxgr is maximum number of grids  many things are
c     dimensioned at, so this is overall. only 1d array
c     though so should suffice. problem is
c     not being able to dimension at maxthreads

      integer locfp_save(maxgr)
      integer locgp_save(maxgr)
      integer lochp_save(maxgr)

c
c  ::::::::::::::; ADVANC :::::::::::::::::::::::::::::::::::::::::::
c  integrate all grids at the input  'level' by one step of its delta(t)
c  this includes:  setting the ghost cells 
c                  advancing the solution on the grid
c                  adjusting fluxes for flux conservation step later
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      mptr = lstart(level)
      hx   = hxposs(level)
      hy   = hyposs(level)
      hz   = hzposs(level)
      delt = possk(level)

c     Lay out linked list into array for easier parallelization
      call prepgrids(listgrids,numgrids(level),level)

c     maxthreads initialized to 1 above in case no openmp
!$    maxthreads = omp_get_max_threads()

c     New code based on 2D
!$OMP PARALLEL DO PRIVATE(j, locnew, locaux, mptr, nx, ny, nz, 
!$OMP&                    mitot, mjtot, mktot,time),
!$OMP&            SHARED(level, nvar, naux, alloc,  delt,
!$OMP&                   nghost, node, rnode, numgrids, listgrids)
!$OMP&            SCHEDULE(dynamic,1),
!$OMP&            DEFAULT(none)
      do j = 1, numgrids(level)
         mptr   = listgrids(j)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         nz     = node(ndkhi,mptr) - node(ndklo,mptr) + 1
         mitot  = nx + 2*nghost
         mjtot  = ny + 2*nghost
         mktot  = nz + 2*nghost
         locnew = node(store1,mptr)
         locaux = node(storeaux,mptr)
         time   = rnode(timemult,mptr)
c     
         call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mktot,
     1              mptr,alloc(locaux),naux)

       end do
!$OMP END PARALLEL DO

c
c save coarse level values if there is a finer level for wave fixup
      if (level+1 .le. mxnest) then
         if (lstart(level+1) .ne. null) then
            call saveqc(level+1,nvar,naux)
         endif
      endif
c
      dtlevnew = rinfinity
      cflmax = 0.d0    !# added by rjl 6/17/05 to keep track of cfl on level

c 
c each grid needs to get extra storage for integration serially (for now at least) in case
c dynamic memory moves the alloc array during expansion.  so get it once and for all for
c biggest chunk it may need. That uses max1d on a side, which for parallel apps.
c is generally much smaller than serial, though this could be a waste.
c remember that max1d includes ghost cells
c
c if necessary can have each thread use max of grids it owns, but then
c cant use dynamic grid assignments.
c
c next loop will get enough storage, it wont be dimensioned as called for
c grids smaller than max1d on a side.
      do j = 1, maxthreads
         
         locfp_save(j) = igetsp(2*max1d*max1d*max1d*nvar)
         locgp_save(j) = igetsp(2*max1d*max1d*max1d*nvar)
         lochp_save(j) = igetsp(2*max1d*max1d*max1d*nvar)

      end do


! Might try using a reduction(min) for dtnew, after basic functionality works
!$OMP PARALLEL DO PRIVATE(j, locold, locnew, mptr,  nx, ny, nz,
!$OMP&                    mitot, mjtot, mktot, locfp, locfm,
!$OMP&                    locgp, locgm, lochp, lochm, ntot, 
!$OMP&                    xlow, ylow, zlow, locaux, lenbc, locsvf,
!$OMP&                    locsvq, locx1d, i, dtnew, time, mythread),
!$OMP&            SHARED(rvol, rvoll, level, nvar, mxnest, alloc, 
!$OMP&                   delt, nghost, intratx, intraty, intratz, hx,
!$OMP&                   hy, hz, naux, listsp, node, rnode, dtlevnew,
!$OMP&                   numgrids,listgrids,
!$OMP&                   locfp_save, locgp_save, lochp_save),
!$OMP&            SCHEDULE (dynamic,1),
!$OMP&            DEFAULT(none)
      do j = 1, numgrids(level)
         mptr   = listgrids(j)
         locold = node(store2, mptr)
         locnew = node(store1, mptr)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         nz     = node(ndkhi,mptr) - node(ndklo,mptr) + 1
         time   = rnode(timemult,mptr)
c
         mitot  = nx + 2*nghost
         mjtot  = ny + 2*nghost

         mktot  = nz + 2*nghost

c        Divide up thread's memory allocation
!$       mythread = omp_get_thread_num()
         locfp = locfp_save(mythread+1)
         locfm = locfp + mitot*mjtot*mktot*nvar
         locgp = locgp_save(mythread+1)
         locgm = locgp + mitot*mjtot*mktot*nvar
         lochp = lochp_save(mythread+1)
         lochm = lochp + mitot*mjtot*mktot*nvar

c
c  copy old soln. values into  next time step's soln. values
c  since integrator will overwrite it. only for grids not at
c  the finest level. finest level grids do not maintain copies
c  of old and new time solution values.
c
         if (level .lt. mxnest) then
            ntot = mitot * mjtot * mktot * nvar
cdir$ ivdep
            do i = 1, ntot
               alloc(locold + i - 1) = alloc(locnew + i - 1)
            end do
         endif
c
         xlow = rnode(cornxlo,mptr) - nghost*hx
         ylow = rnode(cornylo,mptr) - nghost*hy
         zlow = rnode(cornzlo,mptr) - nghost*hz

!$OMP CRITICAL(rv)
         rvol = rvol + nx * ny * nz
         rvoll(level) = rvoll(level) + nx * ny * nz
!$OMP END CRITICAL(rv)


         locaux = node(storeaux,mptr)
c
         if (node(ffluxptr,mptr) .ne. 0) then
            lenbc  = 2*( (nx/intratx(level-1))*(nz/intratz(level-1))
     &                +(ny/intraty(level-1))*(nz/intratz(level-1))
     &                +(nx/intratx(level-1))*(ny/intraty(level-1)))

            locsvf = node(ffluxptr,mptr)
            locsvq = locsvf + nvar*lenbc
            locx1d = locsvq + nvar*lenbc
            call qad(alloc(locnew),mitot,mjtot,mktot,nvar,
     1             alloc(locsvf),alloc(locsvq),lenbc,
     2             intratx(level-1),intraty(level-1),intratz(level-1),
     3             hx,hy,hz,naux,alloc(locaux),alloc(locx1d),delt,mptr)
         endif

c     No gauges in 3D (yet).  Would call dumpgauge here if we had them.

         call stepgrid(alloc(locnew),
     1                alloc(locfm),alloc(locfp),
     2                alloc(locgm),alloc(locgp),
     3                alloc(lochm),alloc(lochp),
     4                mitot,mjtot,mktot,nghost,
     5                delt,dtnew,hx,hy,hz,nvar,
     6                xlow,ylow,zlow,time,mptr,naux,alloc(locaux))

         if (node(cfluxptr,mptr) .ne. 0) then
            call fluxsv(mptr,
     2                  alloc(locfm),alloc(locfp),
     3                  alloc(locgm),alloc(locgp),
     4                  alloc(lochm),alloc(lochp),
     5                  alloc(node(cfluxptr,mptr)),mitot,mjtot,mktot,
     6                  nvar,listsp(level),delt,hx,hy,hz)
         endif
 
         if (node(ffluxptr,mptr) .ne. 0) then
            lenbc = 2*( (nx/intratx(level-1))*(nz/intratz(level-1))
     &                 +(ny/intraty(level-1))*(nz/intratz(level-1))
     &                 +(nx/intratx(level-1))*(ny/intraty(level-1)))
            locsvf = node(ffluxptr,mptr)
            call fluxad(alloc(locfm),alloc(locfp),
     1                  alloc(locgm),alloc(locgp),
     2                  alloc(lochm),alloc(lochp),
     3                  alloc(locsvf),mptr,mitot,mjtot,mktot,nvar,
     4                  lenbc,intratx(level-1),intraty(level-1),
     5                  intratz(level-1),nghost,delt,hx,hy,hz)
         endif

!$OMP CRITICAL (newdt)
         dtlevnew = dmin1(dtlevnew,dtnew)
!$OMP END CRITICAL (newdt)    
c
         rnode(timemult,mptr)  = rnode(timemult,mptr)+delt
      end do
!$OMP END PARALLEL DO

c     new way to reclaim for safety with dynamic memory and openmp
      do j = 1, maxthreads
         call reclam(locfp_save(j),2*max1d*max1d*max1d*nvar)
         call reclam(locgp_save(j),2*max1d*max1d*max1d*nvar)
         call reclam(lochp_save(j),2*max1d*max1d*max1d*nvar)
      end do      

      return
      end
c
c -------------------------------------------------------------
c
       subroutine prepgrids(listgrids, num, level)

       implicit double precision (a-h,o-z)
       include "call.i"
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
