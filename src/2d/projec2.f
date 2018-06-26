c  ::::::::::::::::::::::: PROJEC2 ::::::::::::::::::::::::::::::
! For all newly created fine grids, project area onto a coarser
! grid 2 levels down. Used to recreate grids 1 level down, and
! insure proper level nesting.
! create 
!
! on entry, all coarse grids have already had error estimated, so
! add bad flags.   count number of 'added' flags only.
!
!
!> This subroutine projects all level **level**+2 grids
!! to a level **level** grid and flag the cells being projected as
!! needing refine.
!! In other words, the subroutine modify the flag array of the input
!! grid if part of it is under any grid that is two levels finer.
!!
!! This subroutine is to insure proper level nesting. For example, you
!! just create new level 5 grids. Now you need to ensure that the new
!! level 4 grids encompass the new level 5 grids. To do this, we project
!! level 5 grids to level 3 grids, which means that all locations in 
!! level 3 where level 5 exists are flagged. In this example, level
!! 3 is the **level** parameter on the arguments list of this subroutine.
!! So this subroutine is actually used to ensure **level**+1 grids can
!! encompasses **level**+2 grids.
!!
!! However, note that these cells are flagged with **badpro** parameter
!! defined in [amr_module](@ref amr_module.f90)
!! (not **DOFLAG** as in [flagregions()](@ref flagregions2.f90) and
!! [flag2refine2()](@ref flag2refine2.f90)).
!!
!! **input**: 
!! * **level** 
!! * flag array of the input grid (**rectflags**)
!!
!! **output**: 
!! * **numpro**
!! * flag array of the input grid (**rectflags**)
!!       
!! \param level AMR level of the grid which all fine subgrids are
!! projected onto
!! \param numpro number of additional flagged cells at level **level**
!!               (initialized to 0 in flglvl)
!! \param rectflags array to be flagged 
!! \param ilo global *i* index of the left border of the grid being projected to (being flagged) 
!! \param ihi global *i* index of the right border of the grid being projected to (being flagged) 
!! \param jlo global *j* index of the lower border of the grid being projected to (being flagged) 
!! \param jhi global *i* index of the upper border of the grid being projected to (being flagged) 
!! \param mbuff width of the buffer zone
c  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      subroutine projec2(level,numpro,rectflags,ilo,ihi,jlo,jhi,mbuff)

      use amr_module
      implicit double precision (a-h,o-z)
      dimension rectflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
      logical borderx, bordery
      integer ist(3),iend(3),jst(3),jend(3),ishift(3),jshift(3)

      levpro  =  level + 2
      lrat2x  = intratx(level)*intratx(level+1)
      lrat2y  = intraty(level)*intraty(level+1)

! local variables:
!     mkid    = grid doing the projecting
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
         ! has not intersection
         if (.not.((ixlo .le. ixhi) .and. (jxlo .le. jxhi))) go to 80 
c     
         ! has intersection
         do 60 j = jxlo, jxhi
         do 60 i = ixlo, ixhi
            if (rectflags(i,j) .le. DONTFLAG) then
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
            if (rectflags(ix,jx) .le. DONTFLAG) then
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
                if (rectflags(i,1) .le. DONTFLAG) then
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
                if (rectflags(i,jsize-1) .le. DONTFLAG) then
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
