c
c -------------------------------------------------------------
c
      subroutine ginit(msave, first, nvar, naux, start_time)
c
      use amr_module
      implicit double precision (a-h,o-z)

      logical first
      
      !for setaux timing
      integer :: clock_start, clock_finish, clock_rate
      real(kind=8) :: cpu_start, cpu_finish


c ::::::::::::::::::::::::::::: GINIT ::::::::::::::::::::::::
c
c  initializes soln on all grids at 'level'  by calling qinit
c  if first = true, (first call to init), then allocate the
c  soln storage area too, else was already allocated.
c
c :::::::::::::::::::::::::::::::::::::::;::::::::::::::::::::

      if (msave .eq. 0) go to 99

      level = node(nestlevel,msave)
      hx    = hxposs(level)
      hy    = hyposs(level)
      hz    = hzposs(level)

      mptr  = msave

 10       nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
      nz      = node(ndkhi,mptr) - node(ndklo,mptr) + 1
      mitot   = nx + 2*nghost
      mjtot   = ny + 2*nghost
      mktot  = nz + 2*nghost
      corn1   = rnode(cornxlo,mptr)
      corn2   = rnode(cornylo,mptr)
      corn3   = rnode(cornzlo,mptr)
      if(.not. (first)) go to 20
      loc                 = igetsp(mitot*mjtot*mktot*nvar)
      node(store1,mptr)   = loc
      if (naux .gt. 0) then
         locaux              = igetsp(mitot*mjtot*mktot*naux)
         maxmx = mitot - 2*nghost
         maxmy = mjtot - 2*nghost
         maxmz = mktot - 2*nghost
         mx = maxmx
         my = maxmy
         mz = maxmz
         call system_clocK(clock_start, clock_rate)
         call cpu_time(cpu_start)
         call setaux(nghost,mx,my,mz, corn1,corn2,corn3,
     &               hx,hy,hz,naux,alloc(locaux))
         call system_clock(clock_finish, clock_rate)
         call cpu_time(cpu_finish)

      else
         locaux = 1
      endif
      node(storeaux,mptr) = locaux
      if (level .lt. mxnest) then
         loc2              = igetsp(mitot*mjtot*mktot*nvar)
         node(store2,mptr) = loc2
      endif
      rnode(timemult, mptr) = start_time 
      go to 30
   20 continue
c
c  if 2nd time through, put initial values in store2 so finer grids
c  can be advanced with interpolation of their boundary values.
c  new time soln should still be in location store1.
c
      loc     = node(store2,mptr)
      locaux  = node(storeaux,mptr)
c
   30 continue
c      call qinit(alloc(loc),nvar,mitot,mjtot,mktot,
c     1               corn1-nghost*hx,corn2-nghost*hy,corn3-nghost*hz,
c     2               hx,hy,hz,alloc(locaux),naux)

      maxmx = mitot - 2*nghost
      maxmy = mjtot - 2*nghost
      maxmz = mktot - 2*nghost
      mx = maxmx
      my = maxmy
      mz = maxmz
      call qinit(nvar,nghost,mx,my,mz,corn1,corn2,corn3,
     &           hx,hy,hz,alloc(loc),naux,alloc(locaux))

c
      mptr  = node(levelptr, mptr)
      if (mptr .ne. 0) go to 10
c
c
   99 return
      end
