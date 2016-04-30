c
c ---------------------------------------------------------
c
      subroutine projec2(level,numpro,rectflags,ilo,ihi,jlo,jhi,mbuff)
c
      use amr_module
      implicit double precision (a-h,o-z)
      dimension rectflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
      logical borderx, bordery
      integer ist(3),iend(3),jst(3),jend(3),ishift(3),jshift(3)

c
c  ::::::::::::::::::::::: PROJEC2 ::::::::::::::::::::::::::::::
c  for all newly created fine grids, project area onto a coarser
c  grid 2 levels down. Used to recreate grids 1 level down, and
c  insure proper level nesting.
c
c  on entry, all coarse grids have already had error estimated, so
c  add bad flags.   count number of 'added' flags only.
c
c input parameters:
c    level = project all fine subgrids onto grids at this level.
c output parameters:
c  numpro = number of additional flagged pts. at 'level'.
c           (initialized to 0 in flglvl)
c local variables:
c     mkid    = grid doing the projecting
c  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      levpro  =  level + 2
      lrat2x  = intratx(level)*intratx(level+1)
      lrat2y  = intraty(level)*intraty(level+1)


c
      mkid = newstl(levpro)
 10   if (mkid .eq. 0) go to 90
       ikidlo = node(ndilo,mkid) 
       jkidlo = node(ndjlo,mkid)
       ikidhi = node(ndihi,mkid)
       jkidhi = node(ndjhi,mkid)
c
c  project entire region of fine grids onto rectflag array if intersects
c  possibly take care of buffering.
c  adjust since grid descriptor (integer indices)  is 0 based, 
c  do not projec the buffer region, only interior needs it
c  since buffering will take care of rest (unless ibuff=0-see below)
c
c redo formulas using approach of nestck/baseCheck, simplified to 2 levels
      istc  = ikidlo/intratx(level+1) - 1    ! one level down
      istc  = istc/intratx(level)      - 1    ! project to second level coords
      jstc  = jkidlo/intraty(level+1) - 1
      jstc  = jstc/intraty(level)      - 1
      iendc = ikidhi/intratx(level+1) + 1  
      iendc = iendc/intratx(level)     + 1  
      jendc = jkidhi/intraty(level+1) + 1  
      jendc = jendc/intraty(level)     + 1  

c   if coarse grid not near edge of domain then periodicity wont affect it
      borderx = (istc .le. 0 .or. iendc .ge. iregsz(level)-1)  ! subtract 1 to get last cell index
      bordery = (jstc .le. 0 .or. jendc .ge. jregsz(level)-1)  ! since i/jregsz is num cells

c
c  take care of indices outside actual domain, in non-periodic case first
      if (.not. (xperdom .and. borderx) .and.
     .    .not. (yperdom .and. bordery)) then
         istc  = max(istc,0)
         jstc  = max(jstc,0)
         iendc = min(iendc,iregsz(level))
         jendc = min(jendc,jregsz(level))
         
c  include mbuff in intersection test here since is ok in new alg. to project to buffer region
         ixlo = max(istc, ilo-mbuff)
         ixhi = min(iendc,ihi+mbuff)
         jxlo = max(jstc, jlo-mbuff)
         jxhi = min(jendc,jhi+mbuff)

c        test if coarsened grid mkid intersects with this grids rectflags 
         if (.not.((ixlo .le. ixhi) .and. (jxlo .le. jxhi))) go to 80 
c     
         do 60 j = jxlo, jxhi
         do 60 i = ixlo, ixhi
            if (rectflags(i,j) .eq. goodpt) then
               rectflags(i,j) = badpro
               numpro      = numpro + 1
               if (pprint) write(outunit,101) i,j,mkid
 101           format(' pt.',2i5,' of grid ',i5,' projected' )
            endif
 60      continue
         go to 80            ! done with projected this fine grid in non-periodic case
      endif

c
c periodic case. compute indics on coarsened level to find grids to project to
      call setIndices(ist,iend,jst,jend,iclo,ichi,jclo,jhci,
     .                ishift,jshift,level)

c     compare all regions of coarsened patch with one lbase grid at a time
      do 25 i = 1, 3
         i1 = max(istc,  ist(i))
         i2 = min(iendc, iend(i))
      do 25 j = 1, 3
         j1 = max(jstc,  jst(j))
         j2 = min(jendc, jend(j))

         if (.not. ((i1 .le. i2) .and. (j1 .le. j2))) go to 25
c
c        patch (possibly periodically wrapped) not empty.
c        see if intersects base grid. wrap coords for periodicity
         i1 = i1 + ishift(i)
         i2 = i2 + ishift(i)
         j1 = j1 + jshift(j)
         j2 = j2 + jshift(j)

         ixlo = max(i1,ilo-mbuff)
         ixhi = min(i2,ihi+mbuff)
         jxlo = max(j1,jlo-mbuff)
         jxhi = min(j2,jhi+mbuff)

         if (.not.((ixlo.le.ixhi) .and. (jxlo.le.jxhi))) go to 25

         do jx = jxlo, jxhi
         do ix = ixlo, ixhi
c           project flagged point in intersected regions
            if (rectflags(ix,jx) .eq. goodpt) then
               rectflags(ix,jx) = badpro   ! i,j already coarse grid indices
               numpro      = numpro + 1
               if (pprint) write(outunit,101) ix,jx,mkid
            endif
         end do
         end do

 25   continue
      go to 80   ! down with simple periodic case
c
c repeat above procedure for wrapped area if nec. if ibuff > 0
c this will be caught in shiftset flagging
c  DID NOT MODIFY THIS SPHEREDOM BLOCK WHEN FIXING OTHER BUGS. NEED TO LOOK AT IT
       if (spheredom .and. ibuff .eq. 0) then 
          jstc  = jkidlo/lrat2y
          jendc = jkidhi/lrat2y
          if (jstc .eq. 0) then
             iwrap1 = iregsz(level) - iendc - 1
             iwrap2 = iregsz(level) - istc - 1
c             do 61 i = iwrap1+1, iwrap2+1
             do 61 i = iwrap1, iwrap2  !changing this WITHOUT CHECKING, AS ABOVE. STILL NEED TO CHECK***
                if (rectflags(i,1) .eq. goodpt) then
                  rectflags(i,1) = badpro  ! only need to flag 1 wrapped buffer cell
                  numpro      = numpro + 1
                  if (pprint) write(outunit,101) i,1,mkid
                endif
 61          continue
             
          endif
          if (jendc .eq. jsize-1) then
             iwrap1 = iregsz(level) - iendc - 1
             iwrap2 = iregsz(level) - istc - 1
c             do 62 i = iwrap1+1, iwrap2+1
             do 62 i = iwrap1, iwrap2 !CHANGING W/O CHECKING 
                if (rectflags(i,jsize-1) .eq. goodpt) then
                  rectflags(i,jsize-1) = badpro  ! only need to flag 1 wrapped buffer cell
                  numpro            = numpro + 1
                  if (pprint) write(outunit,101) i,j,mkid
                endif
 62          continue
          endif
       endif
c
c  done with gridpt. loop for grid mkid.
c
 80     mkid = node(levelptr, mkid)
        go to 10
c
 90   if (pprint) then
         write(outunit,102) numpro,level
 102     format(i9,' more pts. projected to level ',i5)

         write(outunit,103) level
 103     format(/,'  from projec: flagged pts. (incl. buffer zone)',
     &           ' at level ',i4,':')

         do 110 j = jhi+mbuff, jlo-mbuff, -1
           write(outunit,104)(int(rectflags(i,j)),i=ilo-mbuff,ihi+mbuff)
104        format(80i1)
 110     continue
      endif
c
 99   return
      end
