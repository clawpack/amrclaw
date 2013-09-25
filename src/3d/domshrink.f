c
c ----------------------------------------------------
c
      subroutine domshrink(iflags2,iflags,idim,jdim,kdim)

      use amr_module
      implicit double precision (a-h, o-z)

      integer*1  iflags2(0:idim+1,0:jdim+1,0:kdim+1)
      integer*1  iflags (0:idim+1,0:jdim+1,0:kdim+1)

c
c :::::::::::::::::::::::::  DOMSHRINK ::::::::::::::::::::::::::::
c
c  shrink domain flags one cell for allowable properly nested domain
c  This is needed even for lcheck = lbase. More shrinking needed
c  for finer levels.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (dprint) then
         write(outunit,*)" from domshrink: on entry, iflags2"
         do 10 kk = 1, kdim
            k = kdim + 1 - kk
            write(outunit,*) 'plane k = ',k
            do 10 jj = 1, jdim
               j = jdim + 1 - jj
               write(outunit,100)(iflags2(i,j,k),i=1,idim)
 100        format(80i1)
 10      continue
      endif

      do 40 k = 1, kdim
      do 40 j = 1, jdim
      do 40 i = 1, idim
      iflags(i,j,k) = iflags2(i,j,k)
      if (iflags2(i  ,j  ,k  ).le.0 .or.
     1    iflags2(i+1,j  ,k  ).le.0 .or. iflags2(i-1,j  ,k  ).le.0 .or.
     2    iflags2(i+1,j+1,k  ).le.0 .or. iflags2(i-1,j+1,k  ).le.0 .or.
     3    iflags2(i  ,j-1,k  ).le.0 .or. iflags2(i  ,j+1,k  ).le.0 .or.
     4    iflags2(i+1,j-1,k  ).le.0 .or. iflags2(i-1,j-1,k  ).le.0 .or.
     .    iflags2(i  ,j  ,k-1).le.0 .or.
     1    iflags2(i+1,j  ,k-1).le.0 .or. iflags2(i-1,j  ,k-1).le.0 .or.
     2    iflags2(i+1,j+1,k-1).le.0 .or. iflags2(i-1,j+1,k-1).le.0 .or.
     3    iflags2(i  ,j-1,k-1).le.0 .or. iflags2(i  ,j+1,k-1).le.0 .or.
     4    iflags2(i+1,j-1,k-1).le.0 .or. iflags2(i-1,j-1,k-1).le.0 .or.
     .    iflags2(i  ,j  ,k+1).le.0 .or.
     1    iflags2(i+1,j  ,k+1).le.0 .or. iflags2(i-1,j  ,k+1).le.0 .or.
     2    iflags2(i+1,j+1,k+1).le.0 .or. iflags2(i-1,j+1,k+1).le.0 .or.
     3    iflags2(i  ,j-1,k+1).le.0 .or. iflags2(i  ,j+1,k+1).le.0 .or.
     4    iflags2(i+1,j-1,k+1).le.0 .or. iflags2(i-1,j-1,k+1).le.0) then
        iflags(i,j,k) = 0
      endif
 40   continue
c
c if border of domain touches a physical boundary then set domain in
c ghost cell as well
c
       if (.not. xperdom) then
         do 55 k = 1, kdim
         do 55 j = 1, jdim
           if (iflags(   1,j,k) .eq. 1) iflags(     0,j,k) = 1
           if (iflags(idim,j,k) .eq. 1) iflags(idim+1,j,k) = 1
 55      continue
       endif
       if (.not. yperdom) then
         do 65 k = 1, kdim
         do 65 i = 1, idim
           if (iflags(i,   1,k) .eq. 1) iflags(i,     0,k) = 1
           if (iflags(i,jdim,k) .eq. 1) iflags(i,jdim+1,k) = 1
 65      continue
       endif
       if (.not. zperdom) then
         do 75 j = 1, jdim
         do 75 i = 1, idim
           if (iflags(i,j,   1) .eq. 1) iflags(i,j,     0) = 1
           if (iflags(i,j,kdim) .eq. 1) iflags(i,j,kdim+1) = 1
 75      continue
       endif

      if (dprint) then
         write(outunit,*)" from domshrink: on exit, iflags"
         do 80 kk = 1, kdim
            k = kdim + 1 - kk
            write(outunit,*) 'plane k = ',k
            do 80 jj = 1, jdim
               j = jdim + 1 - jj
               write(outunit,100)(iflags(i,j,k),i=1,idim)
 80      continue
      endif

      return
      end
