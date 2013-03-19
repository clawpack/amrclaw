c
c ---------------------------------------------------------
c
      subroutine projec2(level,numpro,rectflags,ilo,ihi,jlo,jhi,mbuff)
c
      use amr_module
      implicit double precision (a-h,o-z)
      dimension rectflags(ilo-mbuff:ihi+mbuff,jlo-mbuff:jhi+mbuff)
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
c     iflags  = holds coarser domain flagged points - receives projection
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
c  iflags indexing is 1 based. 
c  do not projec the buffer region, only interior needs it
c  since buffering will take care of rest (unless ibuff=0-see below)
c
      ist  = floor(ikidlo/real(lrat2x))
      jst  = floor(jkidlo/real(lrat2y))
!--      iend = ikidhi/lrat2x
!--      jend = jkidhi/lrat2y
      iend = ceiling((ikidhi+1.d0)/lrat2x) -1
      jend = ceiling((jkidhi+1.d0)/lrat2y) -1

      if (ibuff .eq. 0) then
c     ## ensure proper nesting here, since buffering step won't follow when ibuff 0
        if (ist*lrat2x .eq. ikidlo) ist = ist-1
        if (jst*lrat2y .eq. jkidlo) jst = jst-1
        if ((iend+1)*lrat2x .eq. ikidhi+1) iend = iend+1
        if ((jend+1)*lrat2y .eq. jkidhi+1) jend = jend+1
      endif

      ixlo = max(ist, ilo-mbuff)
      ixhi = min(iend,ihi+mbuff)
      jxlo = max(jst, jlo-mbuff)
      jxhi = min(jend,jhi+mbuff)
      if (.not.((ixlo .le. ixhi) .and. (jxlo .le. jxhi))) go to 80 ! grid mkid doesnt intersect with rectflags
c
c       do 60 j = jst+1, jend+1   !old code, shift indices by 1
c       do 60 i = ist+1, iend+1   ! since iflags used 1-based indexing
c       do 60 j = jst, jend        ! new code into rectflags is 0 based
c       do 60 i = ist, iend       
        do 60 j = jxlo, jxhi
        do 60 i = ixlo, ixhi
           if (rectflags(i,j) .eq. goodpt) then
               rectflags(i,j) = badpro
               numpro      = numpro + 1
               if (pprint) write(outunit,101) i,j,mkid
101            format(' pt.',2i5,' of grid ',i5,' projected' )
           endif
 60    continue
c
c IS THERE SOMETHING TO DO ABOUT PERIODICITY
c
c repeat above procedure for wrapped area if nec. if ibuff > 0
c this will be caught in shiftset flagging
       if (spheredom .and. ibuff .eq. 0) then 
          jst  = jkidlo/lrat2y
          jend = jkidhi/lrat2y
          if (jst .eq. 0) then
             iwrap1 = iregsz(level) - iend - 1
             iwrap2 = iregsz(level) - ist - 1
c             do 61 i = iwrap1+1, iwrap2+1
             do 61 i = iwrap1, iwrap2  !changing this WITHOUT CHECKING, AS ABOVE. STILL NEED TO CHECK***
                if (rectflags(i,1) .eq. goodpt) then
                  rectflags(i,1) = badpro  ! only need to flag 1 wrapped buffer cell
                  numpro      = numpro + 1
                  if (pprint) write(outunit,101) i,1,mkid
                endif
 61          continue
             
          endif
          if (jend .eq. jsize-1) then
             iwrap1 = iregsz(level) - iend - 1
             iwrap2 = iregsz(level) - ist - 1
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
