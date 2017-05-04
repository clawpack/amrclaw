c
c ----------------------------------------------------------------
c
       subroutine setuse(listbc,maxsp,ispot,mkid,
     1                   ilo, ihi,
     2                   iclo,ichi,kflag)
c
c :::::::::::::::::::::::: SETUSE ::::::::::::::::::::::::::::::::
c
c set up boundary list for coarse grid, to be used by fluxsv. 
c loop around boundary of fine grids to do this.  each entry has
c     i, side #, fine grid #, loc in fine grid list for fluxes.
c  for example, side 1 of fine grid fixes side 2 of coarse grid,
c  so coarse grid list will store the # 2.
c  wrt coarse grid, the sides are:
c           1     2       that is, right edge of a coarse cell = 2
c
c  # lkid is the index into the fine grid's saved fluxes.
c  # the fine grid will save all its fluxes all around its
c  # perimeter. lkid tells where the coarse grid should
c  # taking them from. (no ghost cells in this index, but 
c  # it is 1-based for indexing array, not - based for
c  # integer index of grid location).
c
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      use amr_module
      implicit double precision (a-h,o-z)
      parameter(numbcs=4)
      dimension listbc(numbcs,maxsp)


      ibc = ispot
      ist  = iclo - 1
      iend = ichi + 1
c
c  left side (of fine grid, right side of coarse cell)
c
      if (ist .lt. ilo .or. kflag .ne. 1) go to 20
         lkid     = 1

         ispot              = ispot + 1
         listbc(1,ispot)    = ist-ilo+nghost+1
         listbc(2,ispot)    = 2
         listbc(3,ispot)    = mkid
         listbc(4,ispot)    = lkid
         lkid               = lkid + 1
c
c  right side (of fine grid, left of coarse cell)
c (numbered from bottom to top, so not continuous in lkid numbering)
c
 20    if (iend .gt. ihi .or. kflag .ne. 1) go to 30
       lkid     = 2
          ispot              = ispot + 1
          listbc(1,ispot)    = iend-ilo+nghost+1
          listbc(2,ispot)    = 1
          listbc(3,ispot)    = mkid
          listbc(4,ispot)    = lkid
          lkid   = lkid + 1
c
 30    continue
       return
       end
