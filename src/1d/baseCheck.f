c
c ----------------------------------------------------------------
c
       logical function baseCheck(mnew,lbase,ilo,ihi,
     .                            nvar,naux,thisBuff)

       use amr_module
       implicit double precision (a-h, o-z)

       logical debug/.true./
       integer ist(3),iend(3),ishift(3)
       logical borderx
       integer thisBuff

c      index into alloc from iclo:ichi and jclo:jchi, not 0..leni/j. 
       iadd(i) = locm + i - iclo

c ::::::::::::::::::: baseCheck :::::::::::::::::::::::::::
c
c baseCheck - check that potential grid mnew is completely contained
c             in  coarser grids at level 'lbase' (>1) that will
c             stay fixed during this regridding step
c
c   this version tries to do it without using domflags
c   slower but better if cant afford memory over entire domain
c
c   mcheck is one bigger since for proper nesting, cell must be
c   at least one away from boundary of a parent grid, unless
c   on a domain boundary
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::



       levnew = node(nestlevel,mnew)
       borderx = (ilo .eq. 0 .or. ihi .eq. iregsz(levnew)-1)
       

      if (debug) write(outunit,100) mnew,lbase,ilo,ihi,levnew
 100  format("NESTCK1 testing grid ",i5," base level ",i5,/,
     . " new grid from ilo:hi: ",2i12," to ",2i12," at level ",i4)
c
c    on to initializing for the given grid and its nest checking
       levratx = 1
       do 5 lev = lbase, levnew-1
          levratx = levratx * intratx(lev)
 5     continue

c widen by 1 cell (proper nesting), then project to lbase
c this might stick out of domain, fix later
c figure out size for scratch storage on base grid for testing
       iclo = ilo
       ichi = ihi
       do lev = levnew-1,lbase,-1
          iclo = iclo/intratx(lev)
          ichi = ichi/intratx(lev)
          iclo = iclo - 1
          ichi = ichi + 1
          if (debug) then
             write(outunit,111) lev, iclo,ichi
111          format(10x,"at level",i5," projected coords ilo:hi:",2i10,
     .           " jlo:hi:",2i10)
          endif
       end do
c      high end of integer grid index truncates during the divide
c      if it were exactly lined up with coarser grid it would
c      not be properly nested, but since we added one to the index 
c      space, we took care of that already.
       if (debug) then
          write(outunit,108) ilo-1,ihi+1
          write(outunit,109) levratx
 108      format(" enlarged (by 1) fine grid from ilo:hi:",2i12)
 109      format(" refinement factors to base grid of ", 2i12)
          write(outunit,101) iclo,ichi
 101      format("coarsened to lbase, grid from iclo:hi: ",2i12)
       endif

       if (.not. (xperdom .and. borderx)) then
          iclo = max(iclo,0)  ! make sure in domain boundary when checking nesting
          ichi = min(ichi,iregsz(lbase)-1) ! subtract 1 since regsz is number of cells, so -1 is highest index
       endif


       leni = ichi - iclo + 1
       lenrect = leni
       locm = igetsp(lenrect)
       alloc(locm:locm+lenrect-1) = 0.
c
c  if mnew on domain boundary fix flags so ok. 
c  fix extra border, and first/last real edge
       if (ilo .eq. 0 .and. .not. xperdom) then
            alloc(iadd(iclo)) = 1.
            alloc(iadd(iclo+1)) = 1.
       endif
       if (ihi .eq. iregsz(levnew)-1 .and. .not. xperdom) then
             alloc(iadd(ichi)) = 1.
             alloc(iadd(ichi-1)) = 1.
       endif

       mptr = lstart(lbase)
 20       iblo = node(ndilo, mptr) - thisBuff
          ibhi = node(ndihi, mptr) + thisBuff
c
          ! non periodic case, base level coordinates, just mark if nested.
          if (.not. (xperdom .and. borderx)) then
              ixlo = max(iclo,iblo)
              ixhi = min(ichi,ibhi)
              if (.not.(ixlo.le.ixhi)) go to 30
              do ix = ixlo, ixhi
                alloc(iadd(ix))=1.
              end do
              go to 30
          endif
c
c          periodic case:  initialize for potential periodicity
c             each patch divided into 9 regions (some may be empty)
c             e.g. i from (ilo,-1), (0,iregsz(level)-1),(iregsz(level),ihi)
c             except using enlarged grid (ilo-1 to ihi+1)
c
       call setIndices(ist,iend,iclo,ichi,
     .                 ishift,lbase)

c          compare all regions of coarsened patch with one lbase grid at a time
           do 25 i = 1, 3
              i1 = max(iclo,ist(i))
              i2 = min(ichi, iend(i))

              if (.not. (i1 .le. i2)) go to 25
c
c           patch (possibly periodically wrapped) not empty.
c           see if intersects base grid. wrap coords for periodicity
              i1_shifted = i1 + ishift(i)
              i2_shifted = i2 + ishift(i)

              ixlo = max(i1_shifted,iblo)
              ixhi = min(i2_shifted,ibhi)

              if (.not.(ixlo.le.ixhi)) go to 25
c     mark intersected regions with 1
              do ix = ixlo, ixhi
c                need to mark nesting of orig coords, not coarsened shifted indices
                 ix_unshifted = (ix - ishift(i)) ! back to unshifted coords
                                                 ! to mark base grid nesting ok
                 alloc(iadd(ix_unshifted)) = 1.
               end do

 25       continue

 30       mptr = node(levelptr, mptr)
          if (mptr .ne. 0) go to 20

c     output for debugging
      if (debug) then
              write(outunit,344)(int(alloc(iadd(i))), i=iclo,ichi)
 344          format(110i1)
       endif

c
c  if any zeroes left mnew not nested
c
       do 40 i = iclo, ichi
          if (alloc(iadd(i)) .eq. 0) then
             baseCheck = .false.
             go to 99
          endif
 40    continue

c      if made it here then grid is nested
       baseCheck = .true.

 99    call reclam(locm, lenrect)

       return 
       end
