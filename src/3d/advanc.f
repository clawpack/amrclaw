
c --------------------------------------------------------------
c
      subroutine advanc (level,nvar,dtlevnew,vtime,naux)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      logical    vtime

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

 3    continue
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
     1               mptr,alloc(locaux),naux)

	mptr = node(levelptr, mptr)
	if (mptr .ne. 0) go to 3
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
      mptr  = lstart(level)
 5    continue
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
c         ::: get scratch storage for fluxes and slopes
          locfp = igetsp(mitot*mjtot*mktot*nvar)
          locfm = igetsp(mitot*mjtot*mktot*nvar)
          locgp = igetsp(mitot*mjtot*mktot*nvar)
          locgm = igetsp(mitot*mjtot*mktot*nvar)
          lochp = igetsp(mitot*mjtot*mktot*nvar)
          lochm = igetsp(mitot*mjtot*mktot*nvar)
c
c  copy old soln. values into  next time step's soln. values
c  since integrator will overwrite it. only for grids not at
c  the finest level. finest level grids do not maintain copies
c  of old and new time solution values.
c
          if (level .lt. mxnest) then
             ntot   = mitot * mjtot * mktot * nvar
cdir$ ivdep
             do 10 i = 1, ntot
 10            alloc(locold + i - 1) = alloc(locnew + i - 1)
	  endif
c
      xlow = rnode(cornxlo,mptr) - nghost*hx
      ylow = rnode(cornylo,mptr) - nghost*hy
      zlow = rnode(cornzlo,mptr) - nghost*hz
      rvol = rvol + nx * ny * nz
      rvoll(level) = rvoll(level) + nx * ny * nz
      locaux = node(storeaux,mptr)
c
      if (node(ffluxptr,mptr) .ne. 0) then
         lenbc  = 2*( (nx/intratx(level-1))*(nz/intratz(level-1))
     &               +(ny/intraty(level-1))*(nz/intratz(level-1))
     &               +(nx/intratx(level-1))*(ny/intraty(level-1)))
         locsvf = node(ffluxptr,mptr)
         locsvq = locsvf + nvar*lenbc
         locx1d = locsvq + nvar*lenbc
      call qad(alloc(locnew),mitot,mjtot,mktot,nvar,
     1            alloc(locsvf),alloc(locsvq),lenbc,
     2            intratx(level-1),intraty(level-1),intratz(level-1),
     3            hx,hy,hz,naux,alloc(locaux),alloc(locx1d),delt,mptr)
      endif
c
      call stepgrid(alloc(locnew),
     1              alloc(locfm),alloc(locfp),
     2              alloc(locgm),alloc(locgp),
     3              alloc(lochm),alloc(lochp),
     4              mitot,mjtot,mktot,nghost,
     5              delt,dtnew,hx,hy,hz,nvar,
     6              xlow,ylow,zlow,time,mptr,naux,alloc(locaux))

      if (node(cfluxptr,mptr) .ne. 0) then
         call fluxsv(mptr,
     2               alloc(locfm),alloc(locfp),
     3               alloc(locgm),alloc(locgp),
     4               alloc(lochm),alloc(lochp),
     5               alloc(node(cfluxptr,mptr)),mitot,mjtot,mktot,
     6               nvar,listsp(level),delt,hx,hy,hz)
      endif
      if (node(ffluxptr,mptr) .ne. 0) then
         lenbc = 2*( (nx/intratx(level-1))*(nz/intratz(level-1))
     &              +(ny/intraty(level-1))*(nz/intratz(level-1))
     &              +(nx/intratx(level-1))*(ny/intraty(level-1)))
         locsvf = node(ffluxptr,mptr)
         call fluxad(alloc(locfm),alloc(locfp),
     1               alloc(locgm),alloc(locgp),
     2               alloc(lochm),alloc(lochp),
     3               alloc(locsvf),mptr,mitot,mjtot,mktot,nvar,
     4               lenbc,intratx(level-1),intraty(level-1),
     5               intratz(level-1),nghost,delt,hx,hy,hz)
      endif
c
          call reclam(locfp, mitot*mjtot*mktot*nvar)
          call reclam(locfm, mitot*mjtot*mktot*nvar)
          call reclam(locgp, mitot*mjtot*mktot*nvar)
          call reclam(locgm, mitot*mjtot*mktot*nvar)
          call reclam(lochp, mitot*mjtot*mktot*nvar)
          call reclam(lochm, mitot*mjtot*mktot*nvar)
c
          dtlevnew = dmin1(dtlevnew,dtnew)
c
          rnode(timemult,mptr)  = rnode(timemult,mptr)+delt
          mptr            = node(levelptr, mptr)
          if (mptr .ne. 0) go to 5
c
      return
      end
