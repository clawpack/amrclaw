
#define IADD_UP(IVAR,I,J,LOC,NVAR,MITOT) LOC+IVAR-1+NVAR*((J-1)*MITOT+I-1)
#define IADDF_UP(IVAR,I,J,LOCF,NVAR,MI) LOCF+IVAR-1+NVAR*((J-1)*MI+I-1)
#define IADDFAUX_UP(I,J,LOCFAUX,MCAPA,NAUX,MI) LOCFAUX+MCAPA-1+NAUX*((J-1)*MI+(I-1))
#define IADDCAUX_UP(I,J,LOCCAUX,MCAPA,NAUX,MITOT) LOCCAUX+MCAPA-1+NAUX*((J-1)*MITOT+(I-1))

!
! -----------------------------------------------------------
!
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
!
          use amr_module
          implicit double precision (a-h,o-z)


          integer listgrids(numgrids(level))

!$$$  OLD INDEXING
!$$$      iadd(i,j,ivar)  = loc     + i - 1 + mitot*((ivar-1)*mjtot+j-1)
!$$$      iaddf(i,j,ivar) = locf    + i - 1 + mi*((ivar-1)*mj  +j-1)
!$$$      iaddfaux(i,j)   = locfaux + i - 1 + mi*((mcapa-1)*mj + (j-1))
!$$$      iaddcaux(i,j)   = loccaux + i - 1 + mitot*((mcapa-1)*mjtot+(j-1))

!
!
! :::::::::::::::::::::::::: UPDATE :::::::::::::::::::::::::::::::::
! update - update all grids at level 'level'.
!          this routine assumes cell centered variables.
!          the update is done from 1 level finer meshes under it.
! input parameter:
!    level  - ptr to the only level to be updated. levels coarser than
!             this will be at a diffeent time.
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      lget = level
      if (uprint) write(outunit,100) lget
100   format(19h    updating level ,i5)
!     need to set up data structure for parallel distrib of grids
!     call prepgrids(listgrids,numgrids(level),level)

!
!  grid loop for each level
!
      dt     = possk(lget)

!      mptr = lstart(lget)
! 20   if (mptr .eq. 0) go to 85


#ifdef CUDA
!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot, &
!$OMP                     ilo,jlo,ihi,jhi,mkid,iclo,jclo, &
!$OMP                     ichi,jchi,mi,mj,locf,locfaux, &
!$OMP                     iplo,jplo,iphi,jphi,iff,jff,totrat,i,j, &
!$OMP                     ivar,ico,jco,capa,levSt), &
!$OMP          SHARED(lget,numgrids,listgrids,listsp,alloc,nvar,naux, &
!$OMP                    intratx,intraty,nghost,uprint,mcapa,node, &
!$OMP                    listOfGrids,listStart,lstart,level,cflux_hh), &
!$OMP          DEFAULT(none)
#else
!$OMP PARALLEL DO PRIVATE(ng,mptr,loc,loccaux,nx,ny,mitot,mjtot, &
!$OMP                     ilo,jlo,ihi,jhi,mkid,iclo,jclo, &
!$OMP                     ichi,jchi,mi,mj,locf,locfaux, &
!$OMP                     iplo,jplo,iphi,jphi,iff,jff,totrat,i,j, &
!$OMP                     ivar,ico,jco,capa,levSt), &
!$OMP          SHARED(lget,numgrids,listgrids,listsp,alloc,nvar,naux, &
!$OMP                    intratx,intraty,nghost,uprint,mcapa,node, &
!$OMP                    listOfGrids,listStart,lstart,level), &
!$OMP          DEFAULT(none)
#endif
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
!
#ifdef CUDA
         if (associated(cflux_hh(mptr)%ptr) .eq. .false.) go to 25
#else
         if (node(cfluxptr,mptr) .eq. 0) go to 25
#endif

#ifdef CUDA
         call upbnd(cflux_hh(mptr)%ptr,alloc(loc),nvar, &
                    naux,mitot,mjtot,listsp(lget),mptr)
#else
         call upbnd(alloc(node(cfluxptr,mptr)),alloc(loc),nvar, &
                    naux,mitot,mjtot,listsp(lget),mptr)
#endif
!
!  loop through all intersecting fine grids as source updaters.
!
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
!
!  calculate starting and ending indices for coarse grid update, if overlap
!
         iplo = max(ilo,iclo)
         jplo = max(jlo,jclo)
         iphi = min(ihi,ichi)
         jphi = min(jhi,jchi)

         if (iplo .gt. iphi .or. jplo .gt. jphi) go to 75
!
!  calculate starting index for fine grid source pts.
!
         iff    = iplo*intratx(lget) - node(ndilo,mkid) + nghost + 1
         jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
         totrat = intratx(lget) * intraty(lget)
 
         do 71 i = iplo-ilo+nghost+1, iphi-ilo+nghost+1
         do 70 j = jplo-jlo+nghost+1, jphi-jlo+nghost+1
           if (uprint) then
              write(outunit,101) i,j,mptr,iff,jff,mkid
 101          format(' updating pt. ',2i4,' of grid ',i3,' using ',2i4, &
                     ' of grid ',i4)
              write(outunit,102)(alloc(IADD_UP(ivar,i,j,loc,nvar,mitot)), &
                ivar=1,nvar)
 102          format(' old vals: ',4e12.4)
           endif
!
!
!  update using intrat fine points in each direction
!
           do 35 ivar = 1, nvar
 35           alloc(IADD_UP(ivar,i,j,loc,nvar,mitot)) = 0.d0
!
           if (mcapa .eq. 0) then
               do 50 jco  = 1, intraty(lget)
               do 50 ico  = 1, intratx(lget)
               do 40 ivar = 1, nvar
                 alloc(IADD_UP(ivar,i,j,loc,nvar,mitot))= &
                   alloc(IADD_UP(ivar,i,j,loc,nvar,mitot)) + & 
                   alloc(IADDF_UP(ivar,iff+ico-1,jff+jco-1,locf,nvar,mi))
 40              continue
 50            continue
            do 60 ivar = 1, nvar
 60          alloc(IADD_UP(ivar,i,j,loc,nvar,mitot)) = &
                alloc(IADD_UP(ivar,i,j,loc,nvar,mitot))/totrat
               
           else

               do 51 jco  = 1, intraty(lget)
               do 51 ico  = 1, intratx(lget)
               capa = alloc(IADDFAUX_UP(iff+ico-1,jff+jco-1,locfaux,mcapa,naux,mi))
               do 41 ivar = 1, nvar
                 alloc(IADD_UP(ivar,i,j,loc,nvar,mitot)) = &
                    alloc(IADD_UP(ivar,i,j,loc,nvar,mitot)) +  &
                    alloc(IADDF_UP(ivar,iff+ico-1,jff+jco-1,locf,nvar,mi))*capa
 41              continue
 51            continue
            do 61 ivar = 1, nvar
 61          alloc(IADD_UP(ivar,i,j,loc,nvar,mitot)) =  &
                alloc(IADD_UP(ivar,i,j,loc,nvar,mitot))/ &
                (totrat*alloc(IADDCAUX_UP(i,j,loccaux,mcapa,naux,mitot)))
           endif
!
            if (uprint) write(outunit,103)( &
                alloc(IADD_UP(ivar,i,j,loc,nvar,mitot)), ivar=1,nvar)
 103        format(' new vals: ',4e12.4)
!
           jff = jff + intraty(lget)
 70        continue
           iff = iff + intratx(lget)
           jff    = jplo*intraty(lget) - node(ndjlo,mkid) + nghost + 1
 71        continue
!
 75         mkid = node(levelptr,mkid)
            go to 30
!
 80         continue
            end do

!$OMP END PARALLEL DO

!
! 80         mptr = node(levelptr, mptr)
!            go to 20
!
! 85       continue
!
 99   return
      end

