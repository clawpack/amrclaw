c
c ---------------------------------------------------------
c
      subroutine grdfit (lbase,lcheck,nvar,naux,cut,time)
c
      implicit double precision (a-h,o-z)

      include  "call.i"
c
      dimension  corner(nsize,maxcl)
      integer    numptc(maxcl), prvptr
      logical    fit, nestck, cout
      data       cout/.false./
c
c ::::::::::::::::::::: GRDFIT :::::::::::::::::::::::::::::::::;
c  grdfit called by setgrd and regrid to actually fit the new grids
c         on each level. lcheck is the level being error estimated
c         so that lcheck+1 will be the level of the new grids.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      isize = iregsz(lcheck)
      jsize = jregsz(lcheck)
      ksize = kregsz(lcheck)
      ldom2 = igetsp((isize+2)*(jsize+2)*(ksize+2))

c     ## flag all grids at given level based on error ests.
c     ## npts is number of points actually colated - some
c     ## flagged points turned off due to proper nesting requirement.
c     ## (storage based on nptmax calculation however).

      call flglvl (nvar, naux, lcheck, nptmax, index, lbase, ldom2,npts)
      if (npts .eq. 0) go to 99
c
      levnew    = lcheck + 1
      hxfine    = hxposs(levnew)
      hyfine    = hyposs(levnew)
      hzfine    = hzposs(levnew)
c
c     ## call smart_bisect grid gen. to make the clusters
c        till each cluster ok. needs scratch space.
c
       idim = iregsz(lcheck)
       jdim = jregsz(lcheck)
       kdim = kregsz(lcheck)
       lociscr = igetsp(idim+jdim+kdim)
       locjscr = lociscr + idim
       lockscr = locjscr + jdim
       call smartbis(alloc(index),npts,cut,numptc,nclust,lbase,corner,
     1               alloc(lociscr),alloc(locjscr),alloc(lockscr),
     2               idim,jdim,kdim)
       call reclam(lociscr,idim+jdim+kdim)

       if (gprint) then
	  write(outunit,103) nclust
          write(outunit,104) (icl, numptc(icl),icl=1,nclust)
 103      format(' ',i4,' clusters after bisect')
 104      format('         cluster ',i5,' has points: ',i6)
       endif
c
c     ##  for each cluster, fit the actual grid, set up some data structures
c
 50   ibase   =  0
      icl     =  1
      prvptr  =  null
c
 70   mnew      = nodget(dummy)
 75   call  moment(node(1,mnew),alloc(index+numdim*ibase),numptc(icl),
     1             usage)

      if (gprint) write(outunit,100) icl,mnew,usage,numptc(icl)
100   format('         cluster ',i5,' new rect.',i5,
     1       ' usage ',e12.5,' with ',i5,' pts.')

      node(ndilo,mnew)     =  node(ndilo,mnew)*intratx(lcheck)
      node(ndjlo,mnew)     =  node(ndjlo,mnew)*intraty(lcheck)
      node(ndklo,mnew)     =  node(ndklo,mnew)*intratz(lcheck)
      node(ndihi,mnew)     = (node(ndihi,mnew)+1)*intratx(lcheck) - 1
      node(ndjhi,mnew)     = (node(ndjhi,mnew)+1)*intraty(lcheck) - 1
      node(ndkhi,mnew)     = (node(ndkhi,mnew)+1)*intratz(lcheck) - 1
      rnode(cornxlo,mnew)  =  node(ndilo,mnew)*hxfine + xlower
      rnode(cornylo,mnew)  =  node(ndjlo,mnew)*hyfine + ylower
      rnode(cornzlo,mnew)  =  node(ndklo,mnew)*hzfine + zlower
      rnode(cornxhi,mnew)  = (node(ndihi,mnew)+1)*hxfine + xlower
      rnode(cornyhi,mnew)  = (node(ndjhi,mnew)+1)*hyfine + ylower
      rnode(cornzhi,mnew)  = (node(ndkhi,mnew)+1)*hzfine + zlower
      node(nestlevel,mnew) = levnew
      rnode(timemult,mnew) = time
c
c     ##  if new grid doesn't fit in base grid, nestck bisect it
c     ##  and returns 2 clusters where there used to be 1.
c
c 2/28/02 : Add naux to argument list; needed by call to outtre in nestck.
      fit = nestck(mnew,lbase,alloc(index+numdim*ibase),numptc(icl),
     1             numptc,
     2             icl,nclust,alloc(ldom2),isize,jsize,ksize,nvar,naux)
      if (.not. fit) go to 75
c
c     ##  grid accepted. put in list.
      if (newstl(levnew) .eq. null) then
	  newstl(levnew)  = mnew
      else
	  node(levelptr,prvptr) = mnew
      endif
      prvptr = mnew

c     ##  on to next cluster
      ibase  = ibase + numptc(icl)
      icl = icl + 1
      if (icl .le. nclust) go to 70

      if (cout) then
c        2/28/02 : 2d version makes this call to drawrg; What should the
c        3d version do
c        call drawrg(time,lcheck,newstl(levnew),
c        1               nclust,numptc,npts,alloc(index))
      endif


c
c    ##  clean up. for all grids check final size.

      call birect(newstl(levnew))
      call reclam(index, numdim*nptmax)
c
 99   call reclam(ldom2, (isize+2)*(jsize+2)*(ksize+2))

      return
      end
