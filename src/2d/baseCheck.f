c
!> Check that potential grid mnew is completely contained
!! in coarser grids at level **lbase** (>1) that will
!! stay fixed during this regridding step
!! 
!! This version tries to do it without using domflags
!! slower but better if cant afford memory over entire domain
!!
!! For a grid mnew, find the smallest rectangle expressed in 
!! **lbase** index space and compare it to each grid on level **lbase**.
!! If every region in the rectangle has been overlapped during
!! the comparison, the grid mnew is properly nested and the function returns
!! true
!!
!!
c
c ----------------------------------------------------------------
c
       logical function baseCheck(mnew,lbase,ilo,ihi,jlo,jhi,
     .                            nvar,naux,thisBuff)

       use amr_module
       implicit double precision (a-h, o-z)

       logical debug/.false./
       integer ist(3),iend(3),jst(3),jend(3),ishift(3),jshift(3)
       logical borderx, bordery
       integer thisBuff

c      index into alloc from iclo:ichi and jclo:jchi, not 0..leni/j. 
       iadd(i,j) = locm + i - iclo + leni*(j-jclo) 

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
       bordery = (jlo .eq. 0 .or. jhi .eq. jregsz(levnew)-1)
       

c     if (debug) write(outunit,100) mnew,lbase,ilo,ihi,jlo,jhi,levnew
c100  format("NESTCK2 testing grid ",i5," base level ",i5,/,
c    . " new grid from ilo:hi: ",2i12," to ",2i12," at level ",i4)
c
c    on to initializing for the given grid and its nest checking
       levratx = 1
       levraty = 1
       do 5 lev = lbase, levnew-1
          levratx = levratx * intratx(lev)
          levraty = levraty * intraty(lev)
 5     continue

c widen by 1 cell (proper nesting), then project to lbase
c this might stick out of domain, fix later
c figure out size for scratch storage on base grid for testing
       iclo = ilo
       ichi = ihi
       jclo = jlo
       jchi = jhi
       do lev = levnew-1,lbase,-1
          iclo = iclo/intratx(lev)
          ichi = ichi/intratx(lev)
          jclo = jclo/intraty(lev) 
          jchi = jchi/intraty(lev)
          iclo = iclo - 1
          ichi = ichi + 1
          jclo = jclo - 1
          jchi = jchi + 1
c         if (debug) then
c            write(outunit,111) lev, iclo,ichi,jclo,jchi
c111         format(10x,"at level",i5," projected coords ilo:hi:",2i10,
c    .           " jlo:hi:",2i10)
c         endif
       end do
c      high end of integer grid index truncates during the divide
c      if it were exactly lined up with coarser grid it would
c      not be properly nested, but since we added one to the index 
c      space, we took care of that already.
c      if (debug) then
c         write(outunit,108) ilo-1,ihi+1,jlo-1,jhi+1
c         write(outunit,109) levratx,levraty
c108      format(" enlarged (by 1) fine grid from ilo:hi:",2i12,
c    .           " to jlo:hi:", 2i12)
c109      format(" refinement factors to base grid of ", 2i12)
c         write(outunit,101) iclo,ichi,jclo,jchi
c101      format("coarsened to lbase, grid from iclo:hi: ",2i12,
c    .        " to jclo:hi:",2i12)
c      endif

       if (.not. (xperdom .and. borderx) .and. 
     .     .not. (yperdom .and. bordery)) then
          iclo = max(iclo,0)  ! make sure in domain boundary when checking nesting
          jclo = max(jclo,0)  
          ichi = min(ichi,iregsz(lbase)-1) ! subtract 1 since regsz is number of cells, so -1 is highest index
          jchi = min(jchi,jregsz(lbase)-1)
       endif


       leni = ichi - iclo + 1
       lenj = jchi - jclo + 1
       lenrect = leni * lenj
       locm = igetsp(lenrect)
       alloc(locm:locm+lenrect-1) = 0.
c
c  if mnew on domain boundary fix flags so ok. 
c  fix extra border, and first/last real edge
       if (ilo .eq. 0 .and. .not. xperdom) then
         do j = jclo,jchi
            alloc(iadd(iclo,j)) = 1.
            alloc(iadd(iclo+1,j)) = 1.
         end do

       endif
       if (ihi .eq. iregsz(levnew)-1 .and. .not. xperdom) then
          do j = jclo, jchi
             alloc(iadd(ichi,j)) = 1.
             alloc(iadd(ichi-1,j)) = 1.
          end do
       endif
       if (jlo .eq. 0 .and. .not. yperdom) then
         do i = iclo,ichi
            alloc(iadd(i,jclo)) = 1.
            alloc(iadd(i,jclo+1)) = 1.
         end do
       endif
       if (jhi .eq. jregsz(levnew)-1 .and. .not. yperdom) then
          do i = iclo, ichi
             alloc(iadd(i,jchi)) = 1.
             alloc(iadd(i,jchi-1)) = 1.
          end do
       endif

       mptr = lstart(lbase)
 20       iblo = node(ndilo, mptr) - thisBuff
          ibhi = node(ndihi, mptr) + thisBuff
          jblo = node(ndjlo, mptr) - thisbuff
          jbhi = node(ndjhi, mptr) + thisBuff
c
          ! non periodic case, base level coordinates, just mark if nested.
          if ((.not. (xperdom .and. borderx)) .and. 
     .         .not. (yperdom .and. bordery)) then
              ixlo = max(iclo,iblo)
              ixhi = min(ichi,ibhi)
              jxlo = max(jclo,jblo)
              jxhi = min(jchi,jbhi)
              if (.not.((ixlo.le.ixhi) .and. (jxlo.le.jxhi))) go to 30
              do jx = jxlo, jxhi
              do ix = ixlo, ixhi
                alloc(iadd(ix,jx))=1.
              end do
              end do
              go to 30
          endif
c
c          periodic case:  initialize for potential periodicity
c             each patch divided into 9 regions (some may be empty)
c             e.g. i from (ilo,-1), (0,iregsz(level)-1),(iregsz(level),ihi)
c             except using enlarged grid (ilo-1 to ihi+1)
c
       call setIndices(ist,iend,jst,jend,iclo,ichi,jclo,jchi,
     .                 ishift,jshift,lbase)

c          compare all regions of coarsened patch with one lbase grid at a time
           do 25 i = 1, 3
              i1 = max(iclo,ist(i))
              i2 = min(ichi, iend(i))
           do 25 j = 1, 3
              j1 = max(jclo, jst(j))
              j2 = min(jchi, jend(j))

              if (.not. ((i1 .le. i2) .and. (j1 .le. j2))) go to 25
c
c           patch (possibly periodically wrapped) not empty.
c           see if intersects base grid. wrap coords for periodicity
              i1_shifted = i1 + ishift(i)
              i2_shifted = i2 + ishift(i)
              j1_shifted = j1 + jshift(j)
              j2_shifted = j2 + jshift(j)

              ixlo = max(i1_shifted,iblo)
              ixhi = min(i2_shifted,ibhi)
              jxlo = max(j1_shifted,jblo)
              jxhi = min(j2_shifted,jbhi)

              if (.not.((ixlo.le.ixhi) .and. (jxlo.le.jxhi))) go to 25
c     mark intersected regions with 1
              do jx = jxlo, jxhi
              do ix = ixlo, ixhi
c                need to mark nesting of orig coords, not coarsened shifted indices
                 ix_unshifted = (ix - ishift(i)) ! back to unshifted coords
                 jx_unshifted = (jx - jshift(j)) ! to mark base grid nesting ok
                 alloc(iadd(ix_unshifted,jx_unshifted)) = 1.
               end do
               end do

 25       continue

 30       mptr = node(levelptr, mptr)
          if (mptr .ne. 0) go to 20

c     output for debugging
c     if (debug) then
c         do 34 jj = jclo, jchi
c             j = jchi + jclo - jj
c             write(outunit,344)(int(alloc(iadd(i,j))), i=iclo,ichi) 
c344          format(110i1)
c34       continue
c      endif

c
c  if any zeroes left mnew not nested
c
       do 40 j = jclo, jchi
       do 40 i = iclo, ichi
          if (alloc(iadd(i,j)) .eq. 0) then
             baseCheck = .false.
             go to 99
          endif
 40    continue

c      if made it here then grid is nested
       baseCheck = .true.

 99    call reclam(locm, lenrect)

       return 
       end
