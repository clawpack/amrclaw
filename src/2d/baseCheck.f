c
c ----------------------------------------------------------------
c
       logical function baseCheck(mnew,lbase,ilo,ihi,jlo,jhi)

       use amr_module
       implicit double precision (a-h, o-z)

       logical debug/.false./
       integer ist(3),iend(3),jst(3),jend(3),ishift(3),jshift(3)
       ixproj(kk,ll) = (kk + ll*iabs(kk))/ll - iabs(kk)

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



       if (debug) write(outunit,100) mnew,ilo,ihi,jlo,jhi,lbase
 100   format("NESTCK2 testing grid ",i5," base level ",i5,/,
     .  " new grid from ilo:hi: ",2i12," to ",2i12)

       levnew = node(nestlevel,mnew)
c
c initialize for potential periodicity
c each patch divided into 9 regions (some may be empty)
c e.g. i from (ilo,-1), (0,iregsz(level)-1),(iregsz(level),ihi)
c except using enlarged grid (ilo-1 to ihi+1)
c
       call setIndices(ist,iend,jst,jend,ilo-1,ihi+1,jlo-1,jhi+1,
     .                 ishift,jshift,levnew)
c
c    on to initializing for the given grid and its nest checking
       levratx = 1
       levraty = 1
       do 5 lev = lbase, levnew-1
          levratx = levratx * intratx(lev)
          levraty = levraty * intraty(lev)
 5     continue

c widen by 1 cell (proper nesting), then project to lbase
       iclo = ilo-1
       ichi = ihi+1
       jclo = jlo-1
       jchi = jhi+1
       do lev = levnew-1,lbase,-1
          iclo = (dfloat(iclo)/intratx(lev))
          ichi = (dfloat(ichi)/intratx(lev))
          jclo = (dfloat(jclo)/intraty(lev)) 
          jchi = (dfloat(jchi)/intraty(lev))
          if (debug) then
             write(outunit,111) lev, iclo,ichi,jclo,jchi
111          format(10x,"at level",i5," projected coords ilo:hi:",2i10,
     .           " jlo:hi:",2i10)
          endif
       end do
c      high end of integer grid index truncates during the divide
c      if it were exactly lined up with coarser grid it would
c      not be properly nested, but since we added one to the index 
c      space, we took care of that already.
       if (debug) then
          write(outunit,108) ilo-1,ihi+1,jlo-1,jhi+1
          write(outunit,109) levratx,levraty
 108      format(" enlarged fine grid from ilo:hi:",2i12,
     .           " to jlo:hi:", 2i12)
 109      format(" refinement factors to base grid of ", 2i12)
          write(outunit,101) iclo,ichi,jclo,jchi
 101      format("coarsened to lbase, grid from ilo:hi: ",2i12,
     .        " to jlo:hi:",2i12)
       endif

       leni = ichi - iclo + 1
       lenj = jchi - jclo + 1
       lenrect = leni * lenj
       locm = igetsp(lenrect)
c      initialize on projected size of new grid. used to mark nesting         
       do k = 0, lenrect-1   
          alloc(locm + k) = 0.
       end do
c
c      corners dont count for nesting so mark as ok
c  leave 0 for now so match older nestck
!--        alloc(iadd(iclo,jclo)) = 1.
!--        alloc(iadd(iclo,jchi)) = 1.
!--        alloc(iadd(ichi,jclo)) = 1.
!--        alloc(iadd(ichi,jchi)) = 1.
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
 20       iblo = node(ndilo, mptr)
          ibhi = node(ndihi, mptr)
          jblo = node(ndjlo, mptr)
          jbhi = node(ndjhi, mptr)

          if (debug) then
              write(outunit,*)
             write(outunit,102) mptr,lbase,iblo,ibhi,jblo,jbhi
 102        format("TESTING AGAINST GRID ",i5," level ",i5,
     .             " from ilo:hi:",2i10," to jlo:hi:",2i10)
           endif

c          compare all regions of patch with one lbase grid at a time
           do 25 i = 1, 3
              i1 = max(ilo-1,ist(i))
              i2 = min(ihi+1, iend(i))
           do 25 j = 1, 3
              j1 = max(jlo-1, jst(j))
              j2 = min(jhi+1, jend(j))

              if (.not. ((i1 .le. i2) .and. (j1 .le. j2))) go to 25

              if (debug) write(outunit,103)i,j,i1,i2,j1,j2
 103          format("region ",2i10," from ilo:hi:",2i10,
     .               " to jlo:ji:",2i10)
c
c             patch not empty. project coords to level lbase
              iplo = i1+ishift(i)
              iphi = i2+ishift(i)
              jplo = j1+jshift(j)
              jphi = j2+jshift(j)
              do lev = levnew-1,lbase,-1
                iplo = floor(dfloat(iplo)/intratx(lev))
                iphi = ceiling(dfloat(iphi+1)/intratx(lev) - 1)
                jplo = floor(dfloat(jplo)/intraty(lev))
                jphi = ceiling(dfloat(jphi+1)/intraty(lev) - 1)
              end do
              if (debug) then
                 write(outunit,104) i,j,iplo,iphi,jplo,jphi
 104             format("projected coords of region ",2i15," is ",2i12,
     .                  " by ",2i12)
              endif

              ixlo = max(iplo,iblo)
              ixhi = min(iphi,ibhi)
              jxlo = max(jplo,jblo)
              jxhi = min(jphi,jbhi)

              if (.not.((ixlo.le.ixhi) .and. (jxlo.le.jxhi))) go to 25

              if (debug) write(outunit,105) ixlo,ixhi,jxlo,jxhi
 105          format("intersected reg ",2i12," by ",2i12)
c FOR NOW DO NOT HANDLE PERIODIC NEST CHECKING
          
              if (debug) then
                 write(outunit,106) ixlo,ixhi,jxlo,jxhi
                 write(outunit,107) ishift(i),jshift(j)
 106             format("SETTING FROM ",4i5)
 107             format("shifts are ",2i5)
              endif
c     mark intersected regions with 1
              do jx = jxlo, jxhi
              do ix = ixlo, ixhi
c              need to mark nesting of orig coords, not shifted 
               alloc(iadd(ix-ishift(i)/levratx,jx-jshift(j)/levraty))=1.
              end do
              end do

 25       continue

 30       mptr = node(levelptr, mptr)
          if (mptr .ne. 0) go to 20

c     output for debugging
      if (debug) then
          do 34 jj = jclo, jchi
              j = jchi + jclo - jj
              write(outunit,344)(int(alloc(iadd(i,j))), i=iclo,ichi) 
 344          format(110i1)
 34       continue
       endif

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
