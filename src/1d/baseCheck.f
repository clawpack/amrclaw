c
c ----------------------------------------------------------------
c
       logical function baseCheck(mnew,lbase,ilo,ihi)

       use amr_module
       implicit double precision (a-h, o-z)

       logical debug/.false./
       integer ist(3),iend(3),ishift(3)

c      index into alloc from iclo:ichi, not 0..leni.
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



       if (debug) write(outunit,100) mnew,ilo,ihi,lbase
 100   format("NESTCK1 testing grid ",i5," base level ",i5,/,
     .  " new grid from ilo:hi: ",2i12)

       levnew = node(nestlevel,mnew)
c
c initialize for potential periodicity
c each patch divided into 9 regions (some may be empty)
c e.g. i from (ilo,-1), (0,iregsz(level)-1),(iregsz(level),ihi)
c except using enlarged grid (ilo-1 to ihi+1)
c
       call setIndices(ist,iend,ilo-1,ihi+1,
     .                 ishift,levnew)
c
c    on to initializing for the given grid and its nest checking
       levratx = 1
       do 5 lev = lbase, levnew-1
          levratx = levratx * intratx(lev)
 5     continue

c widen by 1 cell (proper nesting), then project to lbase
       iclo = ilo-1
       ichi = ihi+1
       do lev = levnew-1,lbase,-1
          iclo = (dfloat(iclo)/intratx(lev))
          ichi = (dfloat(ichi)/intratx(lev))
          if (debug) then
             write(outunit,111) lev, iclo,ichi
111          format(10x,"at level",i5," projected coords ilo:hi:",2i10)
          endif
       end do
c      high end of integer grid index truncates during the divide
c      if it were exactly lined up with coarser grid it would
c      not be properly nested, but since we added one to the index 
c      space, we took care of that already.
       if (debug) then
          write(outunit,108) ilo-1,ihi+1
          write(outunit,109) levratx
 108      format(" enlarged fine grid from ilo:hi:",2i12)
 109      format(" refinement factors to base grid of ", i12)
          write(outunit,101) iclo,ichi
 101      format("coarsened to lbase, grid from ilo:hi: ",2i12)
       endif

       leni = ichi - iclo + 1
       lenrect = leni
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

       mptr = lstart(lbase)
 20       iblo = node(ndilo, mptr)
          ibhi = node(ndihi, mptr)

          if (debug) then
              write(outunit,*)
             write(outunit,102) mptr,lbase,iblo,ibhi
 102        format("TESTING AGAINST GRID ",i5," level ",i5,
     .             " from ilo:hi:",2i10)
           endif

c          compare all regions of patch with one lbase grid at a time
           do 25 i = 1, 3
              i1 = max(ilo-1,ist(i))
              i2 = min(ihi+1, iend(i))

              if (.not. (i1 .le. i2)) go to 25

              if (debug) write(outunit,103)i,i1,i2
 103          format("region ",i10," from ilo:hi:",2i10)
c
c             patch not empty. project coords to level lbase
              iplo = i1+ishift(i)
              iphi = i2+ishift(i)
              do lev = levnew-1,lbase,-1
                iplo = floor(dfloat(iplo)/intratx(lev))
                iphi = ceiling(dfloat(iphi+1)/intratx(lev) - 1)
              end do
              if (debug) then
                 write(outunit,104) i,iplo,iphi
 104             format("projected coords of region ",i15," is ",2i12)
              endif

              ixlo = max(iplo,iblo)
              ixhi = min(iphi,ibhi)

              if (.not.(ixlo.le.ixhi)) go to 25

              if (debug) write(outunit,105) ixlo,ixhi
 105          format("intersected reg ",2i12)
c FOR NOW DO NOT HANDLE PERIODIC NEST CHECKING
          
              if (debug) then
                 write(outunit,106) ixlo,ixhi
                 write(outunit,107) ishift(i)
 106             format("SETTING FROM ",2i5)
 107             format("shifts are ",i5)
              endif
c     mark intersected regions with 1
              do ix = ixlo, ixhi
c              need to mark nesting of orig coords, not shifted 
               alloc(iadd(ix-ishift(i)/levratx))=1.
              end do

 25       continue

 30       mptr = node(levelptr, mptr)
          if (mptr .ne. 0) go to 20

c     output for debugging
      if (debug) then
          write(outunit,344)(int(alloc(iadd(i))), i=iclo,ichi)
 344      format(110i1)
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
