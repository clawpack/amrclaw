c
c  ----------------------------------------------------------
c
      subroutine domain (nvar,vtime,nx,ny,nz,naux,t0)
c
      implicit double precision (a-h,o-z)
      logical    vtime

      include  "call.i"
c
c  allocate initial coarse grid domain. set node info & initialize grid
c  initial space and time step set here too
c
      mstart = nodget(dummy)
c
c code assumes in many places that lower left corner at (0,0)
c this initial code sets the domain - assumed rectangular
c if it is too large, birect will chop it up into several rectangular
c pieces
c
      rnode(cornxlo,mstart)   = xlower
      rnode(cornylo,mstart)   = ylower
      rnode(cornzlo,mstart)   = zlower
      rnode(cornxhi,mstart)   = xupper
      rnode(cornyhi,mstart)   = yupper
      rnode(cornzhi,mstart)   = zupper
      node(nestlevel,mstart) = 1
      node(levelptr,mstart)  = 0
      lstart(1) = mstart

      if (((nx/2)*2 .ne. nx) .or. (ny/2)*2 .ne. ny
     &                       .or. (nz/2)*2 .ne. nz) then
	 write(outunit,*)" must have even number of cells"
	 write(*,*)      " must have even number of cells"
	 stop
      endif

      node(ndilo,mstart) = 0
      node(ndjlo,mstart) = 0
      node(ndklo,mstart) = 0
      node(ndihi,mstart) = nx-1
      node(ndjhi,mstart) = ny-1
      node(ndkhi,mstart) = nz-1

      lfine = 1
      call  birect(mstart)
      call  ginit (mstart, .true., nvar, naux,t0)
c
c  set stable initial time step using coarse grid data
c
      if (vtime) then
         dt = possk(1)
         dtgrid = dt   !! In case call to estdt doesn't do anything
         mptr = lstart(1)
         dx   = hxposs(1)
         dy   = hyposs(1)
         dz   = hzposs(1)
 60           mitot = node(ndihi,mptr)-node(ndilo,mptr) + 1 + 2*nghost
              mjtot = node(ndjhi,mptr)-node(ndjlo,mptr) + 1 + 2*nghost
              mktot = node(ndkhi,mptr)-node(ndklo,mptr) + 1 + 2*nghost
	      locaux = node(storeaux,mptr)
c 4/1/02 : Added cfl to call to estdt, so that we dont need call.i in estdt
	      call estdt(alloc(node(store1,mptr)),mitot,mjtot,mktot,
     1                   nvar,dx,dy,dz,dtgrid,nghost,alloc(locaux),
     &              naux,cfl)
              dt     = dmin1(dt,dtgrid)
              mptr   = node(levelptr,mptr)
            if (mptr .ne. 0) go to 60
	 possk(1) = dt
      endif
c
c set rest of possk array for refined timesteps
c
      iregsz(1) = nx
      jregsz(1) = ny
      kregsz(1) = nz
      do 70 level = 2, mxnest
	 iregsz(level) = iregsz(level-1) * intratx(level-1)
	 jregsz(level) = jregsz(level-1) * intraty(level-1)
         kregsz(level) = kregsz(level-1) * intratz(level-1)

         possk(level)  = possk(level-1) / dfloat(kratio(level-1))
 70   continue
c
      return
      end
