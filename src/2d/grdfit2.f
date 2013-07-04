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
      logical    fit2, nestck2
      data       cout/.false./
c
c ::::::::::::::::::::: GRDFIT :::::::::::::::::::::::::::::::::;
c  grdfit called by setgrd and regrid to actually fit the new grids
c         on each level. lcheck is the level being error estimated
c         so that lcheck+1 will be the level of the new grids.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c

c ### initialize region start and end indices for new level grids
      iregst(lcheck+1) = iinfinity
      jregst(lcheck+1) = iinfinity
      iregend(lcheck+1) = -1
      jregend(lcheck+1) = -1

c     ## flag all grids at given level based on error ests.
c     ## npts is number of points actually colated - some
c     ## flagged points turned off due to proper nesting requirement.
c     ## (storage based on nptmax calculation however).

      call flglvl2(nvar,naux,lcheck,nptmax,index,lbase,npts,start_time)

      if (npts .eq. 0) go to 99
c
      levnew    = lcheck + 1
      hxfine    = hxposs(levnew)
      hyfine    = hyposs(levnew)
c
c     ## call smart_bisect grid gen. to make the clusters
c        till each cluster ok. needs scratch space.
c
       idim = iregsz(lcheck)
       jdim = jregsz(lcheck)
c       lociscr = igetsp(idim+jdim)
c       locjscr = lociscr + idim
       call smartbis(alloc(index),npts,cut,numptc,nclust,lbase,
     2               corner,idim,jdim)
c     2               corner,alloc(lociscr),alloc(locjscr),idim,jdim)
c       call reclam(lociscr,idim+jdim)

       if (gprint) then
          write(outunit,103) nclust
          write(outunit,104) (icl, numptc(icl),icl=1,nclust)
 103      format(' ',i4,' clusters after bisect')
 104      format('         cluster ',i5,' has points: ',i8)
       endif
c
c     ##  for each cluster, fit the actual grid, set up some data structures
c
 50   ibase   =  0
      icl     =  1
      prvptr  =  null
c
 70   mnew      = nodget(dummy)
c       if (lcheck .eq. 2 .and. (mnew .ne. 6 .and. mnew .ne. 7)) go to 69 
c       if (lcheck .eq. 1 .and. (mnew .ne. 3 .and. mnew .ne. 2 )) go to 69 
 75   call  moment(node(1,mnew),alloc(index+2*ibase),numptc(icl),usage)

      if (gprint) write(outunit,100) icl,mnew,usage,numptc(icl)
100   format('         cluster ',i5,' new rect.',i5,
     1       ' usage ',e12.5,' with ',i5,' pts.')

      node(ndilo,mnew) = node(ndilo,mnew)*intratx(lcheck)
      node(ndjlo,mnew) = node(ndjlo,mnew)*intraty(lcheck)
      node(ndihi,mnew) = (node(ndihi,mnew)+1)*intratx(lcheck) - 1
      node(ndjhi,mnew) = (node(ndjhi,mnew)+1)*intraty(lcheck) - 1
      rnode(cornxlo,mnew)  = node(ndilo,mnew)*hxfine + xlower
      rnode(cornylo,mnew)  = node(ndjlo,mnew)*hyfine + ylower
      rnode(cornxhi,mnew)  = (node(ndihi,mnew)+1)*hxfine + xlower
      rnode(cornyhi,mnew)  = (node(ndjhi,mnew)+1)*hyfine + ylower
      node(nestlevel,mnew)     = levnew
      rnode(timemult,mnew)   = time

      if (gprint) write(outunit,101) (node(i,mnew),i=1,nsize),
     &                              (rnode(i,mnew),i=1,rsize)
 101  format(4i5,4i15,/,4i15,5i15,/,5e15.7)
c
c     ##  if new grid doesn't fit in base grid, nestck bisect it
c     ##  and returns 2 clusters where there used to be 1.
c
c 2/28/02 : Added naux to argument list; needed by call to outtre in nestck

      fit2 = nestck2(mnew,lbase,alloc(index+2*ibase),numptc(icl),numptc,
     1             icl,nclust,nvar, naux)
      if (.not. fit2) go to 75
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
      iregend(levnew) = MAX(iregend(levnew),node(ndihi,mnew))
      jregend(levnew) = MAX(jregend(levnew),node(ndjhi,mnew))

c     ##  on to next cluster
 69     ibase  = ibase + numptc(icl)
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
      if (nptmax .gt. 0) call reclam(index, 2*nptmax) 
c

      return
      end

