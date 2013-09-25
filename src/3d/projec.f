c
c ---------------------------------------------------------
c
      subroutine projec(level,numpro,iflags,isize,jsize,ksize)
c
      use amr_module
      implicit double precision (a-h,o-z)


      integer*1 iflags(0:isize+1,0:jsize+1,0:ksize+1)
c
c  ::::::::::::::::::::::: PROJEC ::::::::::::::::::::::::::::::
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
      lrat2x   = intratx(level)*intratx(level+1)
      lrat2y   = intraty(level)*intraty(level+1)
      lrat2z   = intratz(level)*intratz(level+1)
c
      mkid = newstl(levpro)
 10   if (mkid .eq. 0) go to 90
       ilo = node(ndilo,mkid) 
       jlo = node(ndjlo,mkid)
       klo = node(ndklo,mkid)
       ihi = node(ndihi,mkid)
       jhi = node(ndjhi,mkid)
       khi = node(ndkhi,mkid)
c
c  project entire region of fine grids into iflags array.
c  possibly take care of buffering.
c  adjust since grid descriptor (integer indices)  is 0 based, 
c  iflags indexing is 1 based. 
c
      ist  = ilo/lrat2x 
      jst  = jlo/lrat2y
      kst  = klo/lrat2z
      iend = (ihi+lrat2x-1)/lrat2x 
      jend = (jhi+lrat2y-1)/lrat2y
      kend = (khi+lrat2z-1)/lrat2z
      if (ibuff .eq. 0) then
c     ## ensure proper nesting here, since buffering step won't follow
        if (ist*lrat2x .eq. ilo) ist = ist-1
        if (jst*lrat2y .eq. jlo) jst = jst-1
        if (kst*lrat2z .eq. klo) kst = kst-1
        if (iend*lrat2x .eq. ihi+1) iend = iend+1
        if (jend*lrat2y .eq. jhi+1) jend = jend+1
        if (kend*lrat2z .eq. khi+1) kend = kend+1
      endif
c
       do 60 k = kst+1, kend+1
       do 60 j = jst+1, jend+1
       do 60 i = ist+1, iend+1
           if (iflags(i,j,k) .eq. goodpt) then
               iflags(i,j,k) = badpro
               numpro        = numpro + 1
               if (pprint) write(outunit,101) i,j,k,mkid
101            format(' pt.',3i5,' of grid ',i5,' projected' )
           endif
 60    continue
c
c  done with gridpt. loop for grid mkid.
c
 80     mkid = node(levelptr, mkid)
        go to 10
c
 90   if (numpro .eq. 0) go to 95
      if (pprint) then
        write(outunit,102) numpro,level
 102    format(i7,' more pts. projected to level ',i5)
      endif
c
 95   if (pprint) then
         write(outunit,103) level
 103     format(/,'  from projec: flagged pts. at level ',i4,':')
         do 110 kk = 1, ksize
            k        = ksize + 1 - kk
            write(outunit,*) 'plane k = ',k
         do 110 jj = 1, jsize
            j        = jsize + 1 - jj
            write(outunit,104) (iflags(i,j,k),i=1,isize)
104         format(80i1)
 110     continue
      endif
c
 99   return
      end
