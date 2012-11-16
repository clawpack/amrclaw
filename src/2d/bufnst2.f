c
c -------------------------------------------------------------
c
      subroutine bufnst2 (nvar,naux,numbad,lcheck)
c
      use amr_module
      implicit double precision (a-h,o-z)


      logical    vtime
      integer listgrids(numgrids(lcheck))
      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/
      data       vtime/.false./

c     this indexing is for amrflags array, in flag2refine from 1-mbuff:mx+mbuff
c     but here is from 1:mibuff
      iadd(i,j) = locamrflags + i-1+ mibuff*(j-1)
c
 
c :::::::::::::::::::::::::: BUFNST :::::::::::::::::::::::::::::::::::
c  after error estimation, need to tag the cell for refinement,
c  buffer the tags, take care of level nesting, etc.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
c copy flagged arrays in individual grids (now stored in loctmp)
c into 1 big iflag array (iflags). NO MORE  only do if flag_richardson is true.
c
!$    maxthreads = omp_get_max_threads()
      call prepgrids(listgrids,numgrids(lcheck),lcheck)

      numpro = 0 
      numbad = 0 
      mbuff = max(nghost,ibuff+1)
      time = rnode(timemult,lstart(lcheck))

c      mptr = lstart(lcheck)
c41   continue
!$OMP PARALLEL DO PRIVATE(jg,mptr,ilo,ihi,jlo,jhi,nx,ny,mitot,mjtot),
!$OMP&            PRIVATE(mibuff,mjbuff,locamrflags),
!$OMP&            PRIVATE(loctmp,locbig,j,i,numpro2,numflagged),
!$OMP&            SHARED(numgrids, listgrids,nghost,flag_richardson),
!$OMP&            SHARED(nvar,eprint,maxthreads,mbuff,node),
!$OMP&            SHARED(alloc,lcheck,numpro,mxnest,numbad),
!$OMP&            SCHEDULE (DYNAMIC,1),
!$OMP&            DEFAULT(none)      
      do  jg = 1, numgrids(lcheck)
         mptr = listgrids(jg)
         ilo    = node(ndilo,mptr)
         ihi    = node(ndihi,mptr)
         jlo    = node(ndjlo,mptr)
         jhi    = node(ndjhi,mptr)
         nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
         ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
         mitot  = nx + 2*nghost
         mjtot  = ny + 2*nghost
         mibuff = nx + 2*mbuff
         mjbuff = ny + 2*mbuff
         locamrflags = node(storeflags,mptr)

         if (flag_richardson) then    
            loctmp = node(store2, mptr)
c     call setflags(iflags,isize,jsize,
c     .           alloc(loctmp),nvar,mitot,mjtot,mptr)  
            call addflags(alloc(locamrflags),mibuff,mjbuff,
     .           alloc(loctmp),nvar,mitot,mjtot,mptr)
         endif 
c     still need to reclaim error est space from spest.f 
c     which was saved for possible errest reuse
         locbig = node(tempptr,mptr)
         call reclam(locbig,mitot*mjtot*nvar)
c     
         if (eprint .and. maxthreads .eq. 1) then ! otherwise race for printing
            write(outunit,*)" flagged points before projec2", 
     .           lcheck," grid ",mptr, " (no buff cells)"
            do j = mjbuff-mbuff, mbuff+1, -1
               write(outunit,100)(int(alloc(iadd(i,j))),
     &              i=mbuff+1,mibuff-mbuff)
            enddo
         endif

c     for this version project to each grid separately, no giant iflags
         if (lcheck+2 .le. mxnest) then
            numpro2 = 0
c     write(outunit,*)" calling projec at time ",time
c     write(*,*)" calling projec at time ",time
            call projec2(lcheck,numpro2,alloc(locamrflags),
     .           ilo,ihi,jlo,jhi,mbuff)
c            numpro = numpro + numpro2  not used for now
         endif      

         if (eprint .and. maxthreads .eq. 1) then
            write(outunit,*)" flagged points before buffering on level", 
     .           lcheck," grid ",mptr, " (no buff cells)"
            do 47 j = mjbuff-mbuff, mbuff+1, -1
               write(outunit,100)(int(alloc(iadd(i,j))),
     &              i=mbuff+1,mibuff-mbuff)
 100           format(80i1)
 47         continue
         endif
c     
         if (eprint .and. maxthreads .eq. 1) then
            write(outunit,*)" flagged points after projecting to level", 
     .           lcheck, " grid ",mptr,
     .           "(withOUT buff cells)"
c     .                    "(with buff cells)"
c     buffer zone (wider ghost cell region) now set after buffering
c     so loop over larger span of indices
c     do 49 j = mjbuff, 1, -1
c     write(outunit,100)(int(alloc(iadd(i,j))),i=1,mibuff)
c     49      continue
            do 49 j = mjbuff-mbuff, mbuff+1, -1
               write(outunit,100)(int(alloc(iadd(i,j))),
     .              i=mbuff+1,mibuff-mbuff)
 49         continue
         endif

c
c  diffuse flagged points in all 4 directions to make buffer zones 
c  note that this code flags with a same value as true flagged
c  points, not a different number.
c
c    # first get scratch work space (not that other scratch
c    # arrays above have been reclaimed. 
c
      call shiftset2(alloc(locamrflags),ilo,ihi,jlo,jhi,mbuff)

      if (eprint .and. maxthreads .eq. 1) then
         write(outunit,*)" flagged points after buffering on level", 
     .                    lcheck," grid ",mptr," (WITHOUT buff cells))"
         do 51 j = mjbuff-mbuff, mbuff+1, -1
           write(outunit,100)(int(alloc(iadd(i,j))),
     .                           i=mbuff+1, mibuff-mbuff)
 51      continue
      endif
c   
c   count up
c
      numflagged = 0
      do 82 j = 1, mjbuff
         do 82 i = 1, mibuff
            if (alloc(iadd(i,j)) .ne. goodpt) then 
               numflagged=numflagged + 1
            endif
 82      continue
c     write(outunit,116) numflagged, mptr
 116     format(i5,' points flagged on level ',i4,' grid ',i4)
         node(numflags,mptr) = numflagged
!$OMP CRITICAL(nb)
         numbad = numbad + numflagged   
!$OMP END CRITICAL(nb)

      end do
c     mptr = node(levelptr,mptr)
c     if (mptr .ne. 0) go to 41

      if (verbosity_regrid .gt. 0) then
        write(outunit,*)" total flagged points counted on level ",
     .                  lcheck," is ",numbad
        write(outunit,*)"this may include double counting buffer cells",
     &                  " on  multiple grids"
      endif

      return
      end
