c
c ----------------------------------------------------------
c
      subroutine moment (intrect,badpts,npt,usage)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension     intrect(nsize),badpts(1,npt)
c
c :::::::::::::::::::::::: MOMENT ::::::::::::::::::::::::::::::::::
c  moment = compute enclosing line around flagged points.
c  save some info., even tho. usage might be low and rect. scrapped.

c input parameters:
c     badpts      = x coords of flagged badpts grouped into clusters
c                   are in the first row
c     npt         = num. of badpts. in the cluster.
c
c output parameters:
c     usage       = ratio of flagged to unflagged badpts. in new grid
c                   measures goodness of fit and clustering
c    intrect( )    = stores some info. for grid created herein.
c                   sometimes rect = rnode, sometimes = temp. array.
c                   sometimes intrect = node.
c                   depending on calling prog. (grdfit or expand)
c
c
c :::::::::::::::::::::::: MOMENT ::::::::::::::::::::::::::::::::::
c
      rn = dble(npt)
c
c compute length of enclosing lines to include all flagged badpts.
c
      emx1 = badpts(1,1)
      emn1 = emx1
      do 80 ipt = 1, npt
          if (badpts(1,ipt) .gt. emx1) emx1 = badpts(1,ipt)
          if (badpts(1,ipt) .lt. emn1) emn1 = badpts(1,ipt)
 80   continue
c
c from length of the sides determine corners.
c transform to cell numbers (subtract .5)
c
      intrect(ndilo) = nint(emn1 - .5)
      intrect(ndihi) = nint(emx1 - .5)
c
c compute usage
c
      iside1 = intrect(ndihi) - intrect(ndilo) + 1
      gpall  = iside1
      usage  = rn / gpall
c
      return
      end
