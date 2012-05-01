c
c -----------------------------------------------------------
c
      subroutine colate2 (badpts, len, lcheck, nUniquePts, lbase)
c
      use amr_module
      implicit  double precision (a-h,o-z)
      dimension badpts(2,len)
      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)

c
c    index for flag array now based on integer index space, not 1:mibuff,1:mjbuff
c    but if grid extends outside domain, not supposed to have flagged points
      iadd(i,j) = locamrflags + i-(ilo-mbuff) + mibuff*(j-(jlo-mbuff))
c
c
c *************************************************************
c
c colate2 = takes each grids flagged points at level lcheck
c          and puts their (i,j) cell centered
c          indices into the badpts array.
c          To insure proper nesting, must get rid of flagged points
c          that dont fit into properly nested domain. Grids
c          with flagged points include buffered region (size mbuff)now.
c          THIS NEW VERSION may have duplicate points. need to sort
c          and remove when colating.
c
c if checking flagged pt for nesting is expensive, might consider not doing it
c and revising projec2 instead. if new fine grid not nested, then go through
c flagged points and throw out. But presumably many grids will make it through
c without having to check all points.
c
c *************************************************************
c
c domain flags corresponding to each grid have already been set up.
c colate will check that flagged points nested or throw away
c
       mbuff = max(nghost,ibuff+1)  ! new way of expanding grid to do buffering in place
       index = 0  ! for putting into badpts array


      mptr = lstart(lcheck)
 10      continue
         write(outunit,*)" colating flags on grid ",mptr

c        handle each of 4 sides (in 2D)
c        set tags to negative val. reset to positive if they have a home     
         ilo = node(ndilo,mptr)
         ihi = node(ndihi,mptr)
         jlo = node(ndjlo,mptr)
         jhi = node(ndjhi,mptr)
         nx = ihi - ilo + 1
         ny = jhi - jlo + 1
         mibuff = nx + 2 *mbuff
         mjbuff = ny + 2 *mbuff


         locamrflags = node(storeflags,mptr)
         if (node(numflags,mptr) .eq. 0) go to 70    !simple bypass if no tags

c do we still need setPhysBndry????
         call setPhysBndry(alloc(locamrflags),ilo,ihi,jlo,jhi,
     .                     mbuff,lcheck)
c
         call flagcheck(alloc(locamrflags),ilo,ihi,jlo,jhi,mbuff,
     .                  alloc(node(domflags_upsized,mptr)),mptr)

c
           do 60 j   = jlo-mbuff, jhi+mbuff
           do 60 i   = ilo-mbuff, ihi+mbuff
             if (alloc(iadd(i,j)) .gt. goodpt) then  ! neg means no home was found. throw out
               index = index + 1
c  WARNING: to match orig program note we ADD .5, not subtract. old program used 1 based indexing
c  for grid flagging array. we are using 0 based, so need to add to match previous
c  grid fitting (dont want to change all routines downstream)
c
c PERIODICITY? How tell if wrapped ?
               badpts(1,index) = dble(i)+.5   ! in case periodic, put flagged buffer pt
               badpts(2,index) = dble(j)+.5   ! in badpts in wrapped coords
c              write(outunit,101) badpts(1,index),badpts(2,index)
             else if (alloc(iadd(i,j)) .lt. goodpt) then
               write(outunit,939) i,j
 939         format("NOT NESTED: ignoring point ",2i5)
 101         format(2f6.1)
             endif
 60        continue

 65         continue
 66         continue

c
 70     continue

c  done colating - safe to reclam
        call reclam(locamrflags,mibuff*mjbuff)

        ibytesPerDP = 8
        iflagsize =  (mibuff*mjbuff)/ibytesPerDP+1
        call reclam(node(domflags_base,mptr),iflagsize)
        call reclam(node(domflags_upsized,mptr),iflagsize)

c
        mptr = node(levelptr, mptr)
       if (mptr .ne. 0) go to 10


      npts = index 
      if (gprint) then
        write(outunit,100) npts, lcheck,len
 100    format( i5,' flagged points initially colated on level ',i4,
     .        " badpts len = ",i5)
      endif
c
c colate flagged points into single integer array for quicksorting
c
      call driveSort(npts,badpts,lcheck,nUniquePts,mbuff)
   

 99   return
      end
