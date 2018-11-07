c
c -------------------------------------------------------------
c
      subroutine errest (nvar,naux,lcheck)
c
      use amr_module
      implicit double precision (a-h,o-z)
      
      !for timing setaux
      integer :: clock_start, clock_finish, clock_rate
      real(kind=8) :: cpu_start, cpu_finish

c
c     #   no sense computing new time step if just for error estimation,
c     #   vtime removed from arg list to stepgrid, but new dt not used


c :::::::::::::::::::::::::: ERREST :::::::::::::::::::::::::::::::::::
c for all grids at level lcheck:
c  estimate the error by taking a large (2h,2k) step based on the
c  values in the old storage loc., then take one regular (and for
c  now wasted) step based on the new info.   compare using an
c  error relation for a pth order  accurate integration formula.
c  flag error plane as either bad (needs refinement), or good.
c
c  call the regular integrator on a grid twice as coarse.
c  initialize such a grid directly, instead of trick dimensioning.
c  this is to make other l1 type estimates easier to program.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
       hx   = hxposs(lcheck)
       hy   = hyposs(lcheck)
       hz   = hzposs(lcheck)
       hx2  = 2.d0*hx
       hy2  = 2.d0*hy
       hz2  = 2.d0*hz
       mptr = lstart(lcheck)
 5     continue
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          nz     = node(ndkhi,mptr) - node(ndklo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          mktot  = nz + 2*nghost
          locnew = node(store1,mptr)
          locold = node(store2,mptr)
          locaux = node(storeaux,mptr)
          mi2tot = nx/2  + 2*nghost
          mj2tot = ny/2  + 2*nghost
          mk2tot = nz/2  + 2*nghost
          time   = rnode(timemult,mptr)
          dt     = possk(lcheck)
          tpre   = time - dt
c
c     prepare double the stencil size worth of boundary values,
c            then coarsen them for the giant step integration.
c
          midub = nx+4*nghost
          mjdub = ny+4*nghost
          mkdub = nz+4*nghost
          locdub = igetsp(midub *mjdub *mkdub *(nvar+naux))
          locbgc = igetsp(mi2tot*mj2tot*mk2tot*(nvar+naux))
          node(errptr,mptr) = locbgc
          ngbig = 2*nghost

c         # transfer soln. into grid with twice the ghost cells
          call copysol(alloc(locdub),alloc(locold),nvar,
     &                 mitot,mjtot,mktot,
     1                 nghost,midub,mjdub,mkdub,ngbig)

c
          if (naux .gt. 0) then
              locaxb = locdub + midub *mjdub *mkdub *nvar
              locaxc = locbgc + mi2tot*mj2tot*mk2tot*nvar
              xlowmbc = rnode(cornxlo, mptr)
              ylowmbc = rnode(cornylo, mptr)
              zlowmbc = rnode(cornzlo, mptr)
              xl = xlowmbc - 2*nghost*hx
              yf = ylowmbc - 2*nghost*hy
              zb = zlowmbc - 2*nghost*hz

c	      xl     = rnode(cornxlo, mptr)-2*nghost*hx
c	      yf     = rnode(cornylo, mptr)-2*nghost*hy
c	      zb     = rnode(cornzlo, mptr)-2*nghost*hz
              call system_clock(clock_start, clock_rate)
              call cpu_time(cpu_start)
              call setaux(2*nghost,nx,ny,nz,
     &              xlowmbc,ylowmbc,zlowmbc ,hx,hy,hz,
     &              naux,alloc(locaxb))
              call system_clock(clock_finish, clock_rate)
              call cpu_time(cpu_finish)

              call auxcoarsen(alloc(locaxb),midub ,mjdub ,mkdub ,
     1                        alloc(locaxc),mi2tot,mj2tot,mk2tot,
     2                       naux,auxtype)
          else
              locaxb = 1
          endif

c         # fill it - use enlarged (before coarsening) aux arrays
          call bound(tpre,nvar,
     1               ngbig,alloc(locdub),midub,mjdub,mkdub,mptr,
     2               alloc(locaxb),naux)

c         coarsen by 2 in every direction
          call coarsen(alloc(locdub),midub ,mjdub ,mkdub ,
     1                 alloc(locbgc),mi2tot,mj2tot,mk2tot,nvar)
          call reclam(locdub,midub*mjdub*mkdub*(nvar+naux))
c
c ****************************************************************
c      changed. now locbig filled from spest, always, even if no
c      user estimation routine. reclaim space later tho.
c
c We now fill bndry values at time t = time, in preparation
c for calculating the solution on the grid mptr for error estimation.
c
c          locbig = igetsp(mitot*mjtot*mktot*nvar)
c	        node(tempptr,mptr) = locbig
c         # straight copy into scratch array so don't mess up latest soln.
c          do 10 i = 1, mitot*mjtot*mktot*nvar
c 10          alloc(locbig+i-1) = alloc(locnew+i-1)
c
c          call bound(time,nvar,
c     1               nghost,alloc(locbig),mitot,mjtot,mktot,mptr,
c     1               alloc(locaux),naux)

c ******************************************************************
c
       mptr = node(levelptr,mptr)
       if (mptr .ne. 0) go to 5
c
       hx   = hxposs(lcheck)
       hy   = hyposs(lcheck)
       hz   = hzposs(lcheck)
       hx2  = 2.d0*hx
       hy2  = 2.d0*hy
       hz2  = 2.d0*hz
       dt   = possk(lcheck)
       dt2    = 2. * dt

       mptr = lstart(lcheck)
 25    continue
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          nz     = node(ndkhi,mptr) - node(ndklo,mptr) + 1
          mitot  = nx+ 2*nghost
          mjtot  = ny+ 2*nghost
          mktot  = nz+ 2*nghost
          mi2tot = nx/2 + 2*nghost
          mj2tot = ny/2 + 2*nghost
          mk2tot = nz/2 + 2*nghost
c
c         # this scratch storage will be used both for regular and half
c         # sized grid. different dimensions in stepgrid - do not reuse.
          locfp = igetsp(mitot*mjtot*mktot*nvar)
          locfm = igetsp(mitot*mjtot*mktot*nvar)
          locgp = igetsp(mitot*mjtot*mktot*nvar)
          locgm = igetsp(mitot*mjtot*mktot*nvar)
          lochp = igetsp(mitot*mjtot*mktot*nvar)
          lochm = igetsp(mitot*mjtot*mktot*nvar)
          locaux = node(storeaux,mptr)
c
          locbgc = node(errptr,mptr)
c
          locaxc = locbgc+nvar*mi2tot*mj2tot*mk2tot
c should we set to 1 if naux=0?
c
          evol = evol + (nx/2)*(ny/2)*(nz/2)
          xlow = rnode(cornxlo,mptr) - nghost*hx2
          ylow = rnode(cornylo,mptr) - nghost*hy2
          zlow = rnode(cornzlo,mptr) - nghost*hz2
          call stepgrid(alloc(locbgc),alloc(locfm),alloc(locfp),
     1                alloc(locgm),alloc(locgp),
     &                alloc(lochm),alloc(lochp),
     2                mi2tot,mj2tot,mk2tot,nghost,
     3                dt2,dtnew2,hx2,hy2,hz2,nvar,
     4                xlow,ylow,zlow,tpre,mptr,naux,alloc(locaxc))
c
c  the one giant step based on old values is done. now take
c  one regular step based on new values.
c
          evol   = evol + nx * ny * nz
          xlow   = rnode(cornxlo,mptr) - nghost*hx
          ylow   = rnode(cornylo,mptr) - nghost*hy
          zlow   = rnode(cornzlo,mptr) - nghost*hz
          locbig = node(tempptr,mptr)
          loctmp = node(store2,mptr)
c
c **********************************************************
c *****  changed. spatial error now done from separate routine spest *****
c
c estimate spatial component of error - use old values before
c integration to get accurate boundary gradients
c
c      locerrsp = igetsp(mitot*mjtot*mktot)
c      call errsp(alloc(locbig), alloc(locerrsp), mitot,mjtot,mktot,
c     &           nvar, mptr,nghost, eprint, outunit)
c
c ******************************************************************

      call stepgrid(alloc(locbig),alloc(locfm),alloc(locfp),
     1            alloc(locgm),alloc(locgp),
     &            alloc(lochm),alloc(lochp),
     2            mitot,mjtot,mktot,nghost,
     3            dt,dtnew,hx,hy,hz,nvar,
     4            xlow,ylow,zlow,time,mptr,naux,alloc(locaux))
c
      call reclam(locfp, mitot*mjtot*mktot*nvar)
      call reclam(locfm, mitot*mjtot*mktot*nvar)
      call reclam(locgp, mitot*mjtot*mktot*nvar)
      call reclam(locgm, mitot*mjtot*mktot*nvar)
      call reclam(lochp, mitot*mjtot*mktot*nvar)
      call reclam(lochm, mitot*mjtot*mktot*nvar)
c
      call errf1(alloc(locbig),nvar,
     1          alloc(locbgc),mptr,mi2tot,mj2tot,mk2tot,
     2          mitot,mjtot,mktot,alloc(loctmp))
      call reclam(locbgc  , mi2tot*mj2tot*mk2tot*(nvar+naux))

c      call reclam(locerrsp, mitot *mjtot *mktot)
c reclaiming of locbig now done frombufnst, since this routine will not
c be called if tol<0
c      call reclam(locbig  , mitot *mjtot *mktot *nvar)
c
      mptr = node(levelptr, mptr)
      if (mptr .ne. 0) go to 25

      return
      end
