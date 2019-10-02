c
c -----------------------------------------------------------
c
      subroutine colate1 (badpts, len, lcheck, nUniquePts, lbase)
c
      use amr_module
      implicit  double precision (a-h,o-z)
      dimension badpts(1,len)
      dimension ist(3), iend(3), ishift(3)
      logical db/.false./

c
c    index for flag array now based on integer index space, not 1:mibuff
c    but if grid extends outside domain, not supposed to have flagged points
      iadd(i) = locamrflags + i-(ilo-mbuff)
c
c
c *************************************************************
c
c colate1 = takes each grids flagged points at level lcheck
c          and puts their (i) cell centered
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

c        handle each of 2 sides (in 1D)
c        set tags to negative val. reset to positive if they have a home     
         ilo = node(ndilo,mptr)
         ihi = node(ndihi,mptr)
         nx = ihi - ilo + 1
         mibuff = nx + 2 *mbuff


         locamrflags = node(storeflags,mptr)
         if (node(numflags,mptr) .eq. 0) go to 70    !simple bypass if no tags


c
c  more conservative alg. uses entire buffer in flagging
            imin =  ilo-mbuff
            imax =  ihi+mbuff
            if (.not. xperdom) then
              imin = max(imin,0)
              imax = min(ihi+mbuff,iregsz(lcheck)-1)
            endif
c
c  but to match old alg. use only this one. (not exactly the same? since
c  old alg. used one level?)

c            imin = max(ilo-mbuff,0)
c            imax = min(ihi+mbuff,iregsz(lcheck)-1)

c do we still need setPhysBndry????
c         call setPhysBndry(alloc(locamrflags),ilo,ihi,
c     .                     mbuff,lcheck)
c     pass loop bounds to keep consistent
c     need this next subr. to do integer indexing for iflags
c

         call flagcheck(alloc(locamrflags),ilo,ihi,mbuff,
     .                  alloc(node(domflags2,mptr)),
     .                  imin,imax,mptr)


             do 60 i = imin, imax


c             neg means no home was found. throw out
             if (alloc(iadd(i)) .lt. 0) then
                  write(outunit,939) i
 939              format("NOT NESTED: ignoring point ",2i5)
                  write(*,*)" still have neg points"
                  go to 60
             endif
             if (alloc(iadd(i)) .le. DONTFLAG) go to 60
c
c    got a legit flagged point, bag it.
c
             index = index + 1
c  WARNING: to match orig program note we ADD .5, not subtract. old program used 1 based indexing
c  for grid flagging array. we are using 0 based, so need to add to match previous
c  grid fitting (dont want to change all routines downstream)
c
c  for periodic domains, if flagged pt in buffer zone outside domain
c   wrap it periodically back in before putting on list
              iwrap = i
              if (xperdom) then
                 if (i .lt. 0) iwrap = i + iregsz(lcheck)
                 if (i .ge. iregsz(lcheck)) iwrap = i - iregsz(lcheck)
              endif
c             adding .5 to make it cell centered integer coords
c             note that previous code subtracted .5 since it used 1 based indexing
              badpts(1,index) = dble(iwrap)+.5   ! in case periodic, put flagged buffer pt
             if (db) write(outunit,101) badpts(1,index)
 101          format(2f6.1)

 60        continue

c
 70     continue

c  done colating - safe to reclam
        call reclam(locamrflags,mibuff)

        ibytesPerDP = 8
        iflagsize =  (mibuff)/ibytesPerDP+1
        call reclam(node(domflags_base,mptr),iflagsize)
        call reclam(node(domflags2,mptr),iflagsize)

c
        mptr = node(levelptr, mptr)
       if (mptr .ne. 0) go to 10


      npts = index 
      if (gprint) then
        write(outunit,100) npts, lcheck,len
 100    format( i9,' flagged points initially colated on level ',i4,
     .        " badpts len = ",i10)
      endif
c
c colate flagged points into single integer array for quicksorting
c
c      call drivesort(npts,badpts,lcheck,nUniquePts,mbuff)
c  ### bug found in drivesort- integer overflow
c  ### temp bypass for rnady to run finer grids
      if (lcheck .ge. 6) then
          nUniquePts =    npts
      else
c         # assume points are unique
c         # Cannot check for more than 6 levels due to integer overflow
c         # bug that needs to be fixed!
          call drivesort(npts,badpts,lcheck,nUniquePts,mbuff)
      endif
     

 99   return
      end
