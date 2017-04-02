c
c -----------------------------------------------------------
c
!> Synchronize between all grids on level **level** and grids on
!! level **level**+1.
!! The synchronization includes averaging solution from 
!! level **level**+1 down to level **level** and conservation
!! fix-up near the fine-coarse interface between level **level**
!! grids and level **level**+1 grids.
!!
!! This routine assumes cell centered variables.
!! \param[in] level the only level to be updated (synchronized on). levels coarser than
!! this will be at a diffeent time.
!! \param[in] nvar number of equations for the system
!! \param[in] naux number of auxiliary variables
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
          iadd(ivar,i,j)  = loc    + ivar-1 + nvar*((j-1)*mitot+i-1)
          iaddf(ivar,i,j) = locf   + ivar-1 + nvar*((j-1)*mi+i-1)
          iaddfaux(i,j)   = locfaux + mcapa-1 + naux*((j-1)*mi + (i-1))
          iaddcaux(i,j)   = loccaux + mcapa-1 + naux*((j-1)*mitot+(i-1))
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


!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot,
!$OMP&                    ilo,jlo,ihi,jhi,mkid,iclo,jclo,
!$OMP&                    ichi,jchi,mi,mj,locf,locfaux,
!$OMP&                    iplo,jplo,iphi,jphi,iff,jff,totrat,i,j,
!$OMP&                    ivar,ico,jco,capa,levSt),
!$OMP&         SHARED(lget,numgrids,listgrids,listsp,alloc,nvar,naux,
!$OMP&                   intratx,intraty,nghost,uprint,mcapa,node,
!$OMP&   listOfGrids,listStart,lstart,level),
!$OMP&         DEFAULT(none)

       do ng = 1, numgrids(lget)
         !mptr    = listgrids(ng)
         levSt   = listStart(lget)
         mptr    = listOfGrids(levSt + ng - 1)
         loc     = node(store1,mptr)
         loccaux = node(storeaux,mptr)
         nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mitot   = nx + 2*nghost
         mjtot   = ny + 2*nghost
         ilo     = node(ndilo,mptr)
         jlo     = node(ndjlo,mptr)
         ihi     = node(ndihi,mptr)
         jhi     = node(ndjhi,mptr)
c
         if (node(cfluxptr,mptr) .eq. 0) go to 25
c         locuse = igetsp(mitot*mjtot)
         call upbnd(alloc(node(cfluxptr,mptr)),alloc(loc),nvar,
     1              naux,mitot,mjtot,listsp(lget),mptr)
c     1              mitot,mjtot,listsp(lget),alloc(locuse),mptr)
c         call reclam(locuse,mitot*mjtot)
c
c  loop through all intersecting fine grids as source updaters.
c
 25      mkid = lstart(lget+1)
 30        if (mkid .eq. 0) go to 80
           iclo   = node(ndilo,mkid)/intratx(lget)
           jclo   = node(ndjlo,mkid)/intraty(lget)
           ichi   = node(ndihi,mkid)/intratx(lget)
           jchi   = node(ndjhi,mkid)/intraty(lget)

           mi      = node(ndihi,mkid)-node(ndilo,mkid) + 1 + 2*nghost
           mj      = node(ndjhi,mkid)-node(ndjlo,mkid) + 1 + 2*nghost
           locf    = node(store1,mkid)
           locfaux = node(storeaux,mkid)
c
c  calculate starting and ending indices for coarse grid update, if overlap
c
         iplo = max(ilo,iclo)
         jplo = max(jlo,jclo)
         iphi = min(ihi,ichi)
         jphi = min(jhi,jchi)

         if (iplo .gt. iphi .or. jplo .gt. jphi) go to 75
c
c  calculate starting index for fine grid source pts.
c
         iff    = iplo*intratx(lget) - node(ndilo,mkid) + nghost + 1
         jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
         totrat = intratx(lget) * intraty(lget)
 
         do 71 i = iplo-ilo+nghost+1, iphi-ilo+nghost+1
         do 70 j = jplo-jlo+nghost+1, jphi-jlo+nghost+1
           if (uprint) then
              write(outunit,101) i,j,mptr,iff,jff,mkid
 101          format(' updating pt. ',2i4,' of grid ',i3,' using ',2i4,
     1               ' of grid ',i4)
              write(outunit,102)(alloc(iadd(ivar,i,j)),ivar=1,nvar)
 102          format(' old vals: ',4e12.4)
           endif
c
c
c  update using intrat fine points in each direction
c
           do 35 ivar = 1, nvar
 35           alloc(iadd(ivar,i,j)) = 0.d0
c
           if (mcapa .eq. 0) then
               do 50 jco  = 1, intraty(lget)
               do 50 ico  = 1, intratx(lget)
               do 40 ivar = 1, nvar
                 alloc(iadd(ivar,i,j))= alloc(iadd(ivar,i,j)) + 
     1                        alloc(iaddf(ivar,iff+ico-1,jff+jco-1))
 40              continue
 50            continue
            do 60 ivar = 1, nvar
 60          alloc(iadd(ivar,i,j)) = alloc(iadd(ivar,i,j))/totrat
               
           else

               do 51 jco  = 1, intraty(lget)
               do 51 ico  = 1, intratx(lget)
               capa = alloc(iaddfaux(iff+ico-1,jff+jco-1))
               do 41 ivar = 1, nvar
                 alloc(iadd(ivar,i,j))= alloc(iadd(ivar,i,j)) + 
     1                  alloc(iaddf(ivar,iff+ico-1,jff+jco-1))*capa
 41              continue
 51            continue
            do 61 ivar = 1, nvar
 61          alloc(iadd(ivar,i,j)) = alloc(iadd(ivar,i,j))/
     1                               (totrat*alloc(iaddcaux(i,j)))
           endif
c
            if (uprint) write(outunit,103)(alloc(iadd(ivar,i,j)),
     .                                     ivar=1,nvar)
 103        format(' new vals: ',4e12.4)
c
           jff = jff + intraty(lget)
 70        continue
           iff = iff + intratx(lget)
           jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
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
