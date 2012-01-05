c
c -------------------------------------------------------------
c
      subroutine bufnst2 (nvar,naux,numbad,lcheck,isize,jsize)
c
      use amr_module
      implicit double precision (a-h,o-z)


      logical    vtime
      data       vtime/.false./

      iadd(i,j) = locamrflags + i - ilo-mbuff + mibuff*(j-jlo-mbuff)
c
 
c :::::::::::::::::::::::::: BUFNST :::::::::::::::::::::::::::::::::::
c  after error estimation, need to tag the cell for refinement,
c  buffer the tags, take care of level nesting, etc.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
c copy flagged arrays in individual grids (now stored in loctmp)
c into 1 big iflag array (iflags). NO MORE  only do if tol>0, since otherwise no
c richardson error estimation took place
c
      numpro = 0 
      numbad = 0 
      mbuff = max(nghost,ibuff)

      mptr = lstart(lcheck)
 41   continue
      ilo    = node(ndilo,mptr)
      ihi    = node(ndihi,mptr)
      jlo    = node(ndjlo,mptr)
      jhi    = node(ndjhi,mptr)
      nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
      mitot  = nx + 2*nghost
      mjtot  = ny + 2*nghost
      locamrflags = node(storeflags,mptr)
c     negative tol means richardson not done      
      if (tol .gt. 0) then    
            loctmp = node(store2, mptr)
c            call setflags(iflags,isize,jsize,
c     .           alloc(loctmp),nvar,mitot,mjtot,mptr)  
            mibuff = nx + 2*mbuff
            mjbuff = ny + 2*mbuff
            call addflags(alloc(locamrflags),mibuff,mjbuff,
     .                    alloc(loctmp),nvar,mitot,mjtot,mptr)
      endif 
c     still need to reclaim error est space from spest.f 
c     which was saved for possible errest reuse
      locbig = node(tempptr,mptr)
      call reclam(locbig,mitot*mjtot*nvar)
c
c for this version project to each grid separately, no giant iflags
      if (lcheck+2 .le. mxnest) then
         call projec2(lcheck,numpro2,alloc(locamrflags),
     .                ilo,ihi,jlo,jhi,mbuff)
         numpro = numpro + numpro2
      endif      

       if (eprint) then
         write(outunit,*)" flagged points before buffering on level", 
     .                   lcheck," grid ",mptr, "(no buff cells)"
         do 47 j = jhi, jlo, -1
           write(outunit,100)(alloc(iadd(i,j)),i=ilo,ihi)
 100       format(80i1)
 47      continue
      endif
c
c  project finer grids to insure level nesting
c      numpro = 0
c      if (lcheck+2 .le. mxnest) then
c         call projec(lcheck,numpro,iflags,isize,jsize)
c      endif

      if (eprint) then
         write(outunit,*)" flagged points after projecting to level", 
     .                    lcheck, " grid ",mptr
         write(outunit,*) " with ",numpro," additional points projected"
         do 49 j = jhi+mbuff, jlo-mbuff, -1
           write(outunit,100)(alloc(iadd(i,j)),i=ilo-mbuff,ihi+mbuff)
 49      continue
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

      if (eprint) then
         write(outunit,*)" flagged points after buffering on level", 
     .                    lcheck," grid ",mptr
         do 48 j = jhi+mbuff, jlo-mbuff, -1
           write(outunit,100)(alloc(iadd(i,j)),i=ilo-mbuff, ihi+mbuff)
 48      continue
      endif
c   
c   count up
c
          numflagged = 0
          do 82 j = jlo-mbuff, jhi+mbuff
          do 82 i = ilo-mbuff, ihi+mbuff
           if (alloc(iadd(i,j)) .ne. goodpt) numflagged=numflagged + 1
 82       continue
          write(outunit,116) numflagged, mptr
 116      format(i5,' points flagged on level ',i4,' grid ',i4)
          node(numflags,mptr) = numflagged
          numbad = numbad + numflagged

      mptr = node(levelptr,mptr)
      if (mptr .ne. 0) go to 41

      write(outunit,*)" total points flagged on level ",lcheck,
     .                " is ",numbad
      write(outunit,*)"this may include double counting of buffer cells"

      return
      end
