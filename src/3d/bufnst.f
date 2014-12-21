c
c -------------------------------------------------------------
c
      subroutine bufnst (nvar,naux,numbad,lcheck,iflags,
     .                   isize,jsize,ksize)
c
      use amr_module
      implicit double precision (a-h,o-z)


      integer*1  iflags (0:isize+1,0:jsize+1,0:ksize+1)
      logical    vtime
      data       vtime/.false./
c
 
c :::::::::::::::::::::::::: BUFNST :::::::::::::::::::::::::::::::::::
c  after error estimation, need to tag the cell for refinement,
c  buffer the tags, take care of level nesting, etc.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
c copy flagged arrays in individual grids (now stored in loctmp)
c into 1 big iflag array (iflags). only do if tol>0, since otherwise no
c richardson error estimation took place
c
      mptr = lstart(lcheck)
 41   continue
      nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
      nz     = node(ndkhi,mptr) - node(ndklo,mptr) + 1
      mitot  = nx + 2*nghost
      mjtot  = ny + 2*nghost
      mktot  = nz + 2*nghost
c     negative tol means richardson not done      
      if (tol .gt. 0) then    
         loctmp = node(store2, mptr)
         call setflags(iflags,isize,jsize,ksize,
     .                 alloc(loctmp),nvar,mitot,mjtot,mktot,mptr)
      endif 
c     still need to reclaim error est space from spest.f which was saved for possible errest reuse
      locbig = node(tempptr,mptr)
      call reclam(locbig,mitot*mjtot*mktot*nvar)
      mptr = node(levelptr,mptr)
      if (mptr .ne. 0) go to 41


      if (eprint) then
         write(outunit,*)" flagged points before buffering on level", 
     .                   lcheck
         do 47 kk = 1, ksize
           k = ksize + 1 - kk
         do 47 jj = 1, jsize
           j = jsize + 1 - jj
           write(outunit,100)(iflags(i,j,k),i=1,isize)
 100       format(80i1)
 47      continue
      endif
c
c  project finer grids to insure level nesting
      numpro = 0
      if (lcheck+2 .le. mxnest) then
         call projec(lcheck,numpro,iflags,isize,jsize,ksize)
      endif

      if (eprint) then
         write(outunit,*)" flagged points after projecting to level", 
     .                    lcheck
         write(outunit,*) " with ",numpro," additional points projected"
         do 49 kk = 1, ksize
           k = ksize + 1 - kk
         do 49 jj = 1, jsize
           j = jsize + 1 - jj
           write(outunit,100)(iflags(i,j,k),i=1,isize)
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
c      ldom3 = igetsp((isize+2)*(jsize+2)*(ksize+2))
c
!--      do 55 inum = 1, ibuff
!--
!--          call shiftset(iflags,alloc(ldom3),+1, 0, 0,isize,jsize,ksize)
!--          call shiftset(iflags,alloc(ldom3),-1, 0, 0,isize,jsize,ksize)
!--          call shiftset(iflags,alloc(ldom3), 0,+1, 0,isize,jsize,ksize)
!--          call shiftset(iflags,alloc(ldom3), 0,-1, 0,isize,jsize,ksize)
!--          call shiftset(iflags,alloc(ldom3), 0, 0,+1,isize,jsize,ksize)
!--          call shiftset(iflags,alloc(ldom3), 0, 0,-1,isize,jsize,ksize)
!--
!-- 55   continue
      call shiftset(iflags, isize,jsize,ksize)

      if (eprint) then
         write(outunit,*)" flagged points after buffering on level", 
     .                    lcheck
         do 48 kk = 1, ksize
           k = ksize + 1 - kk
         do 48 jj = 1, jsize
           j = jsize + 1 - jj
           write(outunit,100)(iflags(i,j,k),i=1,isize)
 48      continue
      endif
c   
c   count up
c
       numbad = 0 
       do 82 i = 1, isize
       do 82 j = 1, jsize
       do 82 k = 1, ksize
         if (iflags(i,j,k) .ne. goodpt) numbad = numbad + 1
 82    continue
       if (eprint) write(outunit,116) numbad, lcheck
 116   format(i9,' points flagged on level ',i4)

c      call reclam(ldom3,(isize+2)*(jsize+2)*(ksize+2))

      return
      end
