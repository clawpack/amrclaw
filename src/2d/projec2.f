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
       ilo = node(ndilo,mkid) 
       jlo = node(ndjlo,mkid)
       ihi = node(ndihi,mkid)
       jhi = node(ndjhi,mkid)
c
c  project entire region of fine grids onto rectflag array if intersects
c  possibly take care of buffering.
c  adjust since grid descriptor (integer indices)  is 0 based, 
c  iflags indexing is 1 based. 
c
      ist  = ilo/lrat2x 
      jst  = jlo/lrat2y
      iend = ihi/lrat2x
      jend = jhi/lrat2y
      if (ibuff .eq. 0) then
c     ## ensure proper nesting here, since buffering step won't follow
        if (ist*lrat2x .eq. ilo) ist = ist-1
        if (jst*lrat2y .eq. jlo) jst = jst-1
        if ((iend+1)*lrat2x .eq. ihi+1) iend = iend+1
        if ((jend+1)*lrat2y .eq. jhi+1) jend = jend+1
      endif
c
       do 60 j = jst+1, jend+1
       do 60 i = ist+1, iend+1
           if (rectflags(i,j) .eq. goodpt) then
               rectflags(i,j) = badpro
               numpro      = numpro + 1
               if (pprint) write(outunit,101) i,j,mkid
101            format(' pt.',2i5,' of grid ',i5,' projected' )
           endif
 60    continue
c
c repeat above procedure for wrapped area if nec. if ibuff > 0
c this will be caught in shiftset flagging
       if (spheredom .and. ibuff .eq. 0) then 
          jst  = jlo/lrat2y
          jend = jhi/lrat2y
          if (jst .eq. 0) then
             iwrap1 = iregsz(level) - iend - 1
             iwrap2 = iregsz(level) - ist - 1
             do 61 i = iwrap1+1, iwrap2+1
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
             do 62 i = iwrap1+1, iwrap2+1
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
 90   if (numpro .eq. 0) go to 95
      write(outunit,102) numpro,level
 102  format(i7,' more pts. projected to level ',i5)
c
 95   if (pprint) then
         write(outunit,103) level
 103     format(/,'  from projec: flagged pts. at level ',i4,':')
         do 110 jj = 1, mjbuff
            j        = mjbuff + 1 - jj
            write(outunit,104) (rectflags(i,j),i=ilo-mbuff,iihi+mbuff)
104         format(80i1)
 110     continue
      endif
c
 99   return
      end
