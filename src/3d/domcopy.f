c
c ----------------------------------------------------
c
      subroutine domcopy(iflags2,iflags,isize,jsize,ksize)

      use amr_module
      implicit double precision (a-h, o-z)

      integer*1  iflags2(0:isize+1,0:jsize+1,0:ksize+1)
      integer*1  iflags (0:isize+1,0:jsize+1,0:ksize+1)

c
c ::::::::::::::::::::::::::: DOMCOPY :::::::::::::::::::::
c 
c  domain flags are in iflags. copy into iflags2.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::


      do 10 k = 0, ksize+1
      do 10 j = 0, jsize+1
      do 10 i = 0, isize+1
        iflags2(i,j,k) = iflags(i,j,k)
 10   continue
c
c  take care of periodicity again
c
      if (xperdom) then
         do 35 k = 0, ksize+1
         do 35 j = 0, jsize+1
           iflags2(      0,j,k) = iflags2(isize,j,k)
           iflags2(isize+1,j,k) = iflags2(    1,j,k)
 35      continue
       else
         do 65 k = 1, ksize
         do 65 j = 1, jsize
           if (iflags2(    1,j,k) .eq. 1) iflags2(      0,j,k) = 1
           if (iflags2(isize,j,k) .eq. 1) iflags2(isize+1,j,k) = 1
 65      continue
      endif
      if (yperdom) then
         do 45 k = 0, ksize+1
         do 45 i = 0, isize+1
           iflags2(i,      0,k) = iflags2(i,jsize,k)
           iflags2(i,jsize+1,k) = iflags2(i,    1,k)
 45      continue
       else
         do 55 k = 1, ksize
         do 55 i = 1, isize
           if (iflags2(i,    1,k) .eq. 1) iflags2(i,      0,k) = 1
           if (iflags2(i,jsize,k) .eq. 1) iflags2(i,jsize+1,k) = 1
 55      continue
      endif
      if (zperdom) then
         do 75 j = 0, jsize+1
         do 75 i = 0, isize+1
           iflags2(i,j,      0) = iflags2(i,j,ksize)
           iflags2(i,j,ksize+1) = iflags2(i,j,    1)
 75      continue
       else
         do 85 j = 1, jsize
         do 85 i = 1, isize
           if (iflags2(i,j,    1) .eq. 1) iflags2(i,j,      0) = 1
           if (iflags2(i,j,ksize) .eq. 1) iflags2(i,j,ksize+1) = 1
 85      continue
       end if

      if (dprint) then
         write(outunit,*)" from domcopy: domflags "
         do 40 kk = 1, ksize
         k = ksize + 1 - kk
         write(outunit,*) 'plane k = ',k
         do 40 jj = 1, jsize
         j = jsize + 1 - jj
         write(outunit,100)(iflags2(i,j,k),i=1,isize)
 100     format(80i1)
 40      continue
      endif

 
      return
      end
