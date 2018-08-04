c
!> Takes flagged points on all grids on level lcheck
!! and pack their (i,j) cell centered
!! indices into the badpts array.
!! Points in the badpts array are unique and sorted based on
!! one dimensional packing of their 2D indices.
c
c -----------------------------------------------------------
c
      subroutine colate2 (badpts, len, lcheck, nUniquePts, lbase)
c
      use amr_module
      implicit  double precision (a-h,o-z)
      dimension badpts(2,len)
      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)
      logical db/.false./
      integer*8 largestIntEquiv,ifac1,ifac2
      integer   largestIntEquiv_default

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
c        write(outunit,*)" colating flags on grid ",mptr

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


c
c  more conservative alg. uses entire buffer in flagging
            jmin =  jlo-mbuff
            jmax =  jhi+mbuff
            imin =  ilo-mbuff
            imax =  ihi+mbuff
            if (.not. xperdom) then
              imin = max(imin,0)
              imax = min(ihi+mbuff,iregsz(lcheck)-1)
            endif
            if (.not. yperdom) then
              jmin = max(jmin,0)
              jmax = min(jhi+mbuff,jregsz(lcheck)-1)
            endif
c
c  but to match old alg. use only this one. (not exactly the same? since
c  old alg. used one level?)

c            jmin = max(jlo-mbuff,0)
c            jmax = min(jhi+mbuff,jregsz(lcheck)-1)
c            imin = max(ilo-mbuff,0)
c            imax = min(ihi+mbuff,iregsz(lcheck)-1)

c do we still need setPhysBndry????
c         call setPhysBndry(alloc(locamrflags),ilo,ihi,jlo,jhi,
c     .                     mbuff,lcheck)
c     pass loop bounds to keep consistent
c     need this next subr. to do integer indexing for iflags
c
         call flagcheck(alloc(locamrflags),ilo,ihi,jlo,jhi,mbuff,
     .                  alloc(node(domflags2,mptr)),
     .                  imin,imax,jmin,jmax,mptr)


             do 60 j = jmin, jmax
             do 60 i = imin, imax

c             neg means no home was found. throw out
             if (alloc(iadd(i,j)) .lt. 0) then
                  write(outunit,939) i,j
 939              format("NOT NESTED: ignoring point ",2i5)
                  write(*,*)" still have neg points"
                  go to 60
             endif
             if (alloc(iadd(i,j)) .le. DONTFLAG) go to 60
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
              jwrap = j
              if (yperdom) then
                 if (j .lt. 0) jwrap = j + jregsz(lcheck)
                 if (j .ge. jregsz(lcheck)) jwrap = j - jregsz(lcheck)
              endif
c             adding .5 to make it cell centered integer coords
c             note that previous code subtracted .5 since it used 1 based indexing
              badpts(1,index) = dble(iwrap)+.5   ! in case periodic, put flagged buffer pt
              badpts(2,index) = dble(jwrap)+.5   ! in badpts in wrapped coords
             if (db) write(outunit,101) badpts(1,index), badpts(2,index)
 101          format(2f6.1)

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
c     sorting uses one dimensional packing of 2D indices
c     check if domain will fit  in integer*4.  (largestSingle approx 2**30)
c     if not, just leave the duplicate, but rememer that efficiency
c     of grids won't be correct (since divide by number of flaged pts in grid)
c     If necessary, do whole process in integer*8 - then will have enough
c     room, but will have to convert quicksort routine and drivesort
c     the variable largestIntEquiv already declared integer*8 above.
      ifac1 = iregsz(lcheck)
      ifac2 = jregsz(lcheck)
      largestIntEquiv =  ifac1+mbuff +
     .             (ifac1+2*mbuff)*(ifac2+mbuff)
      largestIntEquiv_default =  iregsz(lcheck)+mbuff +
     .             (iregsz(lcheck)+2*mbuff)*(jregsz(lcheck)+mbuff)

c     ! if get different answer with extra precision then bypass sorting alg.
      if (largestIntEquiv .ne. largestIntEquiv_default) then
c       ## sorting alg will have integer overflow
c       ## just use all flagged points in making grids
c       ## this means "efficiency" count will be incorrect for
c       ## this and higher levels
          nUniquePts =  npts  ! bad name - they are not unique
      else
          call drivesort(npts,badpts,lcheck,nUniquePts,mbuff)
      endif


 99   return
      end
