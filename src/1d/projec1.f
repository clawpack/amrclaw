c
c ---------------------------------------------------------
c
      subroutine projec1(level,numpro,rectflags,ilo,ihi,mbuff)
c
      use amr_module
      implicit double precision (a-h,o-z)
      dimension rectflags(ilo-mbuff:ihi+mbuff)
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
c
      mkid = newstl(levpro)
 10   if (mkid .eq. 0) go to 90
       ikidlo = node(ndilo,mkid)
       ikidhi = node(ndihi,mkid)
c
c  project entire region of fine grids onto rectflag array if intersects
c  possibly take care of buffering.
c  adjust since grid descriptor (integer indices)  is 0 based, 
c  iflags indexing is 1 based. 
c  do not projec the buffer region, only interior needs it
c  since buffering will take care of rest (unless ibuff=0-see below)
c
      ist  = floor(ikidlo/real(lrat2x))
!--      iend = ikidhi/lrat2x
      iend = ceiling((ikidhi+1.d0)/lrat2x) -1

      if (ibuff .eq. 0) then
c     ## ensure proper nesting here, since buffering step won't follow when ibuff 0
        if (ist*lrat2x .eq. ikidlo) ist = ist-1
        if ((iend+1)*lrat2x .eq. ikidhi+1) iend = iend+1
      endif

      ixlo = max(ist, ilo-mbuff)
      ixhi = min(iend,ihi+mbuff)
      if (.not.(ixlo .le. ixhi)) go to 80 ! grid mkid doesnt intersect with rectflags
c
c       !old code, shift indices by 1
c       do 60 i = ist+1, iend+1   ! since iflags used 1-based indexing
c       do 60 i = ist, iend        ! new code into rectflags is 0 based
        do 60 i = ixlo, ixhi
           if (rectflags(i) .le. DONTFLAG) then
               rectflags(i) = badpro
               numpro      = numpro + 1
               if (pprint) write(outunit,101) i,mkid
101            format(' pt.',2i5,' of grid ',i5,' projected' )
           endif
 60    continue
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

           write(outunit,104)(int(rectflags(i)),i=ilo-mbuff,ihi+mbuff)
104        format(80i1)
      endif
c
 99   return
      end
