c
c ------------------------------------------------------------------
c
       subroutine setuse(listbc,maxsp,ispot,mkid,
     1                   ilo,ihi,jlo,jhi,klo,khi,
     2                   iclo,ichi,jclo,jchi,kclo,kchi,nghost)
c
c ::::::::::::::::::::::::: SETUSE ::::::::::::::::::::::::::::::::
c
c set up boundary list for coarse grid, to be used by fluxsv. 
c loop around boundary of c fine grids to do this.  each entry has
c     i, j, k, side #, fine grid #, loc in fine grid list for fluxes.
c  for example, side 1 of fine grid fixes side 3 of coarse grid,
c  so coarse grid list will store the # 3.
c  wrt coarse grid, in the xy plane, the sides are:
c              2
c           1     3            that is, right edge of a coarse cell = 3
c              4                         rear edge of a coarse cell = 2
c  and in the xz plane the sides are:
c              6
c           1     3            that is, right face of a coarse cell = 3
c              5                       bottom face of a coarse cell = 5
c
c  # lkid is the index into the fine grid's saved fluxes.
c  # the fine grid will save all its fluxes all around its
c  # perimeter. lkid tells where the coarse grid should
c  # taking them from.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      implicit double precision (a-h,o-z)
      parameter(numbcs=6)
      dimension listbc(numbcs,maxsp)

      ibc = ispot
      ist  = iclo - 1
      iend = ichi + 1
      jst  = jclo - 1
      jend = jchi + 1
      kst  = kclo - 1
      kend = kchi + 1
c
c  left side
c
       if (ist .lt. ilo) go to 20
c2D    lkid     = max(jlo,jclo) - jclo + 1
       do 15 k  = max(klo,kclo), min(khi,kchi)
       lkid     = (k-kclo)*(jchi-jclo+1) +
     1            max(jlo,jclo) - jclo + 1
       do 10 j  = max(jlo,jclo), min(jhi,jchi)
          ispot              = ispot + 1
          listbc(1,ispot)    = ist-ilo+nghost+1
          listbc(2,ispot)    = j-jlo+nghost+1
          listbc(3,ispot)    = k-klo+nghost+1
          listbc(4,ispot)    = 3
          listbc(5,ispot)    = mkid
          listbc(6,ispot)    = lkid
          lkid               = lkid + 1
 10    continue
 15    continue
c
c   rear side
c
 20    if (jend .gt. jhi) go to 40
c2D    lkid       = (jchi-jclo+1) + max(ilo,iclo)-iclo + 1
       do 35 k  = max(klo,kclo), min(khi,kchi)
       lkid     = (kchi-kclo+1)*(jchi-jclo+1) +
     1            (k-kclo)*(ichi-iclo+1) +
     2            max(ilo,iclo)-iclo + 1
       do 30 i    = max(ilo,iclo), min(ihi,ichi)
          ispot              = ispot + 1
          listbc(1,ispot)    = i-ilo+nghost+1
          listbc(2,ispot)    = jend-jlo+nghost+1
          listbc(3,ispot)    = k-klo+nghost+1
          listbc(4,ispot)    = 4
          listbc(5,ispot)    = mkid
          listbc(6,ispot)    = lkid
          lkid               = lkid + 1
 30    continue
 35    continue
c
c  right side (numbered from bottom to top, so not continuous)
c
 40    if (iend .gt. ihi) go to 60
c2D    lkid     = (ichi-iclo+1)+(jchi-jclo+1)
c2D  .               + max(jlo,jclo) - jclo + 1
       do 55 k  = max(klo,kclo), min(khi,kchi)
       lkid     = (kchi-kclo+1)*(jchi-jclo+1) +
     1            (kchi-kclo+1)*(ichi-iclo+1) +
     2            (k-kclo)*(jchi-jclo+1) +
     3            max(jlo,jclo) - jclo + 1
       do 50 j  = max(jlo,jclo), min(jhi,jchi)
          ispot              = ispot + 1
          listbc(1,ispot)    = iend-ilo+nghost+1
          listbc(2,ispot)    = j-jlo+nghost+1
          listbc(3,ispot)    = k-klo+nghost+1
          listbc(4,ispot)    = 1
          listbc(5,ispot)    = mkid
          listbc(6,ispot)    = lkid
          lkid   = lkid + 1
 50    continue
 55    continue
c
c  front  side (numbered left to right, so not continuous)
c
 60    if (jst .lt. jlo) go to 80
c2D    lkid   =  2*(jchi-jclo+1) + (ichi-iclo+1) + max(ilo,iclo)-iclo+1
       do 75 k  = max(klo,kclo), min(khi,kchi)
       lkid     = (kchi-kclo+1)*(jchi-jclo+1) +
     1            (kchi-kclo+1)*(ichi-iclo+1) +
     2            (kchi-kclo+1)*(jchi-jclo+1) +
     3            (k-kclo)*(ichi-iclo+1) +
     4            max(ilo,iclo)-iclo + 1
       do 70 i  = max(ilo,iclo), min(ihi,ichi)
          ispot              = ispot + 1
          listbc(1,ispot)    = i-ilo+nghost+1
          listbc(2,ispot)    = jst-jlo+nghost+1
          listbc(3,ispot)    = k-klo+nghost+1
          listbc(4,ispot)    = 2
          listbc(5,ispot)    = mkid
          listbc(6,ispot)    = lkid
          lkid   = lkid + 1
 70    continue
 75    continue
c
c  bottom side
c
 80    if (kst .lt. klo) go to 100
       do 95 j  = max(jlo,jclo), min(jhi,jchi)
       lkid     = (kchi-kclo+1)*(jchi-jclo+1) +
     1            (kchi-kclo+1)*(ichi-iclo+1) +
     2            (kchi-kclo+1)*(jchi-jclo+1) +
     3            (kchi-kclo+1)*(ichi-iclo+1) +
     4            (j-jclo)*(ichi-iclo+1) +
     5            max(ilo,iclo)-iclo + 1
       do 90 i  = max(ilo,iclo), min(ihi,ichi)
          ispot              = ispot + 1
          listbc(1,ispot)    = i-ilo+nghost+1
          listbc(2,ispot)    = j-jlo+nghost+1
          listbc(3,ispot)    = kst-klo+nghost+1
          listbc(4,ispot)    = 6
          listbc(5,ispot)    = mkid
          listbc(6,ispot)    = lkid
          lkid   = lkid + 1
 90    continue
 95    continue
c
c   top side
c
100    if (kend .gt. khi) go to 120
       do 115 j = max(jlo,jclo), min(jhi,jchi)
       lkid     = (kchi-kclo+1)*(jchi-jclo+1) +
     1            (kchi-kclo+1)*(ichi-iclo+1) +
     2            (kchi-kclo+1)*(jchi-jclo+1) +
     3            (kchi-kclo+1)*(ichi-iclo+1) +
     4            (jchi-jclo+1)*(ichi-iclo+1) +
     5            (j-jclo)*(ichi-iclo+1) +
     6            max(ilo,iclo)-iclo + 1
       do 110 i = max(ilo,iclo), min(ihi,ichi)
          ispot              = ispot + 1
          listbc(1,ispot)    = i-ilo+nghost+1
          listbc(2,ispot)    = j-jlo+nghost+1
          listbc(3,ispot)    = kend-klo+nghost+1
          listbc(4,ispot)    = 5
          listbc(5,ispot)    = mkid
          listbc(6,ispot)    = lkid
          lkid               = lkid + 1
 110   continue
 115   continue
c
 120   continue
       return
       end
