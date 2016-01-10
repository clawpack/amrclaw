c
c ---------------------------------------------------------
c
      subroutine grdfit (lbase,lcheck,nvar,naux,cut,time,start_time)
c
      use amr_module
      implicit double precision (a-h,o-z)

c
      dimension  corner(nsize,maxcl)
      integer    numptc(maxcl), prvptr
      logical    fit, nestck, cout
      data       cout/.false./
      integer*1  i1flags(iregsz(lcheck)+2,jregsz(lcheck)+2,
     .                                    kregsz(lcheck)+2)
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
c      ldom2 = igetsp((isize+2)*(jsize+2)*(ksize+2))

c ### initialize region start and end indices for new level grids
      iregst(lcheck+1) = iinfinity
      jregst(lcheck+1) = iinfinity
      kregst(lcheck+1) = iinfinity
      iregend(lcheck+1) = -1
      jregend(lcheck+1) = -1
      kregend(lcheck+1) = -1

c     ## flag all grids at given level based on error ests.
c     ## npts is number of points actually colated - some
c     ## flagged points turned off due to proper nesting requirement.
c     ## (storage based on nptmax calculation however).

      call flglvl (nvar, naux, lcheck, nptmax, index, lbase, 
     .             i1flags,npts,start_time,isize,jsize,ksize)
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
       call smartbis(alloc(index),npts,cut,numptc,nclust,lbase,corner,
     2               idim,jdim,kdim)

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
 70   mnew      = nodget()
 75   call  moment(node(1,mnew),alloc(index+3*ibase),numptc(icl),
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
      fit = nestck(mnew,lbase,alloc(index+3*ibase),numptc(icl),
     1            numptc,icl,nclust,i1flags,isize,jsize,ksize,nvar,naux)
      if (.not. fit) go to 75
c
c     ##  grid accepted. put in list.
      if (newstl(levnew) .eq. null) then
         newstl(levnew)  = mnew
      else
         node(levelptr,prvptr) = mnew
      endif
      prvptr = mnew
c     # keep track of min and max location of grids at this level
      iregst(levnew)  = MIN(iregst(levnew), node(ndilo,mnew))
      jregst(levnew)  = MIN(jregst(levnew), node(ndjlo,mnew))
      kregst(levnew)  = MIN(kregst(levnew), node(ndklo,mnew))
      iregend(levnew) = MAX(iregend(levnew),node(ndihi,mnew))
      jregend(levnew) = MAX(jregend(levnew),node(ndjhi,mnew))
      kregend(levnew) = MAX(kregend(levnew),node(ndkhi,mnew))

c     ##  on to next cluster
      ibase  = ibase + numptc(icl)
      icl = icl + 1
      if (icl .le. nclust) go to 70

c
c    ##  clean up. for all grids check final size.

      call birect(newstl(levnew))
 99   continue

c    ## may have npts 0 but array was allocated due to initially flagged points
c    ## that were not allowed for proper nesting or other reasons. in this case
c    ## the array was still allocated, so need to test further to see if colating
c    ## array space needs to be reclaimed
      if (nptmax .gt. 0)  call reclam(index, 3*nptmax)
c
c 99   call reclam(ldom2, (isize+2)*(jsize+2)*(ksize+2))

      return
      end
