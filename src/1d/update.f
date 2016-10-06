c
c -----------------------------------------------------------
c
      subroutine update (level, nvar, naux)
c
      use amr_module
      implicit double precision (a-h,o-z)


      integer listgrids(numgrids(level))

c$$$  OLD INDEXING
c$$$      iadd(i,j,ivar)  = loc     + i - 1 + mitot*((ivar-1)*mjtot+j-1)
c$$$      iaddf(i,j,ivar) = locf    + i - 1 + mi*((ivar-1)*mj  +j-1)
c$$$      iaddfaux(i,j)   = locfaux + i - 1 + mi*((mcapa-1)*mj + (j-1))
c$$$      iaddcaux(i,j)   = loccaux + i - 1 + mitot*((mcapa-1)*mjtot+(j-1))

c   NEW INDEXING, ORDER SWITCHED
      iadd(ivar,i)  = loc    + ivar-1 + nvar*(i-1)
      iaddf(ivar,i) = locf   + ivar-1 + nvar*(i-1)
      iaddfaux(i)   = locfaux + mcapa-1 + naux*(i-1)
      iaddcaux(i)   = loccaux + mcapa-1 + naux*(i-1)
c
c
c :::::::::::::::::::::::::: UPDATE :::::::::::::::::::::::::::::::::
c update - update all grids at level 'level'.
c          this routine assumes cell centered variables.
c          the update is done from 1 level finer meshes under it.
c input parameter:
c    level  - ptr to the only level to be updated. levels coarser than
c             this will be at a diffeent time.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      lget = level
      if (uprint) write(outunit,100) lget
100   format(19h    updating level ,i5)
c     need to set up data structure for parallel distrib of grids
c     call prepgrids(listgrids,numgrids(level),level)

c
c  grid loop for each level
c
      dt     = possk(lget)

c      mptr = lstart(lget)
c 20   if (mptr .eq. 0) go to 85


!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,mitot,
!$OMP&                    ilo,ihi,mkid,iclo,
!$OMP&                    ichi,mi,locf,locfaux,
!$OMP&                    iplo,iphi,iff,totrat,i,
!$OMP&                    ivar,ico,capa,levSt),
!$OMP&         SHARED(lget,numgrids,listgrids,listsp,alloc,nvar,naux,
!$OMP&                   intratx,nghost,uprint,mcapa,node,
!$OMP&   listOfGrids,listStart,lstart,level),
!$OMP&         DEFAULT(none)

       do ng = 1, numgrids(lget)
         !mptr    = listgrids(ng)
         levSt   = listStart(lget)
         mptr    = listOfGrids(levSt + ng - 1)
         loc     = node(store1,mptr)
         loccaux = node(storeaux,mptr)
         nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
         mitot   = nx + 2*nghost
         ilo     = node(ndilo,mptr)
         ihi     = node(ndihi,mptr)
c
         if (node(cfluxptr,mptr) .eq. 0) go to 25
c         locuse = igetsp(mitot*mjtot)
         call upbnd(alloc(node(cfluxptr,mptr)),alloc(loc),nvar,
     1              naux,mitot,listsp(lget),mptr)
c     1              mitot,mjtot,listsp(lget),alloc(locuse),mptr)
c         call reclam(locuse,mitot*mjtot)
c
c  loop through all intersecting fine grids as source updaters.
c
 25      mkid = lstart(lget+1)
 30        if (mkid .eq. 0) go to 80
           iclo   = node(ndilo,mkid)/intratx(lget)
           ichi   = node(ndihi,mkid)/intratx(lget)

           mi      = node(ndihi,mkid)-node(ndilo,mkid) + 1 + 2*nghost
           locf    = node(store1,mkid)
           locfaux = node(storeaux,mkid)
c
c  calculate starting and ending indices for coarse grid update, if overlap
c
         iplo = max(ilo,iclo)
         iphi = min(ihi,ichi)

         if (iplo .gt. iphi) go to 75
c
c  calculate starting index for fine grid source pts.
c
         iff    = iplo*intratx(lget) - node(ndilo,mkid) + nghost + 1
         totrat = intratx(lget)
 
         do 71 i = iplo-ilo+nghost+1, iphi-ilo+nghost+1
           if (uprint) then
              write(outunit,101) i,mptr,iff,mkid
 101          format(' updating pt. ',i4,' of grid ',i3,' using ',i4,
     1               ' of grid ',i4)
              write(outunit,102)(alloc(iadd(ivar,i)),ivar=1,nvar)
 102          format(' old vals: ',4e12.4)
           endif
c
c
c  update using intrat fine points in each direction
c
           do 35 ivar = 1, nvar
 35           alloc(iadd(ivar,i)) = 0.d0
c
           if (mcapa .eq. 0) then
               do 50 ico  = 1, intratx(lget)
               do 40 ivar = 1, nvar
                 alloc(iadd(ivar,i))= alloc(iadd(ivar,i)) +
     1                        alloc(iaddf(ivar,iff+ico-1))
 40              continue
 50            continue
            do 60 ivar = 1, nvar
 60          alloc(iadd(ivar,i)) = alloc(iadd(ivar,i))/totrat
               
           else

               do 51 ico  = 1, intratx(lget)
               capa = alloc(iaddfaux(iff+ico-1))
               do 41 ivar = 1, nvar
                 alloc(iadd(ivar,i))= alloc(iadd(ivar,i)) +
     1                  alloc(iaddf(ivar,iff+ico-1))*capa
 41              continue
 51            continue
            do 61 ivar = 1, nvar
 61          alloc(iadd(ivar,i)) = alloc(iadd(ivar,i))/
     1                               (totrat*alloc(iaddcaux(i)))
           endif
c
            if (uprint) write(outunit,103)(alloc(iadd(ivar,i)),
     .                                     ivar=1,nvar)
 103        format(' new vals: ',4e12.4)
c

           iff = iff + intratx(lget)
 71        continue
c
 75         mkid = node(levelptr,mkid)
            go to 30
c
 80         continue
            end do

!$OMP END PARALLEL DO

c
c 80         mptr = node(levelptr, mptr)
c            go to 20
c
c 85       continue
c
 99   return
      end
