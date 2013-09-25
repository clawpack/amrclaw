c
c ----------------------------------------------------
c
      subroutine domup(iflags2,iflags,ibase,jbase,kbase,
     1                                isize,jsize,ksize,lev)

      use amr_module
      implicit double precision (a-h, o-z)

      integer*1  iflags2(0:isize+1,0:jsize+1,0:ksize+1)
      integer*1  iflags (0:ibase+1,0:jbase+1,0:kbase+1)

c
c ::::::::::::::::::::::::::: DOMUP :::::::::::::::::::::
c 
c  domain flags are in iflags. copy into iflags2, allowing
c  for change of level and dimension
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (dprint) then
         write(outunit,*)" from domup: domflags (before expansion)"
         do 5 kk = 1, kbase
         k = kbase + 1 - kk
         write(outunit,*) 'plane k = ',k
         do 5 jj = 1, jbase
         j = jbase + 1 - jj
         write(outunit,100)(iflags(i,j,k),i=1,ibase)
 5       continue
      endif

      do 10 k = 0, ksize+1
      do 10 j = 0, jsize+1
      do 10 i = 0, isize+1
         iflags2(i,j,k) = 0
 10   continue

      do 20 k = 1, kbase
         kfine = (k-1) * intratz(lev)
      do 20 j = 1, jbase
         jfine = (j-1) * intraty(lev)
      do 20 i = 1, ibase
         ifine = (i-1) * intratx(lev)
         do 25 mk = 1, intratz(lev)
         do 25 mj = 1, intraty(lev)
         do 25 mi = 1, intratx(lev)
            iflags2(ifine+mi,jfine+mj,kfine+mk) = iflags(i,j,k)  
 25      continue
 20   continue
c
c  take care of periodicity again or if border of domain touches a 
c  physical boundary then set domain in ghost cell as well
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
 65    continue
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
      endif
c
c     corner stuff
c
      do k = 1, ksize
        if (iflags2(0,1,k)+iflags2(1,0,k).eq.2)  iflags2(0,0,k) = 1
        if (iflags2(0,jsize,k)+iflags2(1,jsize+1,k).eq.2)
     .                                       iflags2(0,jsize+1,k) = 1
        if (iflags2(isize,0,k)+iflags2(isize+1,1,k).eq.2)
     .                                       iflags2(isize+1,0,k) = 1
        if (iflags2(isize,jsize+1,k)+iflags2(isize+1,jsize,k).eq.2)
     .                                 iflags2(isize+1,jsize+1,k) = 1
      end do
      do i = 1, isize
        if (iflags2(i,0,1)+iflags2(i,1,0).eq.2)  iflags2(i,0,0) = 1
        if (iflags2(i,jsize,0)+iflags2(i,jsize+1,0).eq.2)
     .                                       iflags2(i,jsize+1,0) = 1
        if (iflags2(i,0,ksize)+iflags2(i,1,ksize+1).eq.2)
     .                                       iflags2(i,0,ksize+1) = 1
        if (iflags2(i,jsize+1,ksize)+iflags2(i,jsize,ksize+1).eq.2)
     .                                 iflags2(i,jsize+1,ksize+1) = 1
      end do
      do j = 1, jsize
        if (iflags2(0,j,1)+iflags2(1,j,0).eq.2)  iflags2(0,j,0) = 1
        if (iflags2(isize,j,0)+iflags2(isize+1,j,1).eq.2)
     .                                       iflags2(isize+1,j,0) = 1
        if (iflags2(0,j,ksize)+iflags2(1,j,ksize+1).eq.2)
     .                                       iflags2(0,j,ksize+1) = 1
        if (iflags2(isize+1,j,ksize)+iflags2(isize,j,ksize+1).eq.2)
     .                                 iflags2(isize+1,j,ksize+1) = 1
      end do

c
c     The neigbors of a corner are the three locations each having an
c     address which differs from that of the corner by 1 in just one of
c     the three indexes. If a corner index is 0, the corresponding neighbor
c     index may be 0 or 1. If a corner index is size+1, the corresponding
c     neighbor index may be size+1 or size (where size => [i,j,k]size).
c
      if (( iflags2(      0,      0,      1)
     &     +iflags2(      0,      1,      0)
     &     +iflags2(      1,      0,      0)) .eq. 3)
     .      iflags2(      0,      0,      0)   =   1
 
      if (( iflags2(isize+1,      0,      1)
     &     +iflags2(isize+1,      1,      0)
     &     +iflags2(isize  ,      0,      0)) .eq. 3)
     .      iflags2(isize+1,      0,      0)   =   1
 
      if (( iflags2(      0,jsize+1,      1)
     &     +iflags2(      0,jsize  ,      0)
     &     +iflags2(      1,jsize+1,      0)) .eq. 3)
     .      iflags2(      0,jsize+1,      0)   =   1
 
      if (( iflags2(isize+1,jsize+1,      1)
     &     +iflags2(isize+1,jsize  ,      0)
     &     +iflags2(isize  ,jsize+1,      0)) .eq. 3)
     .      iflags2(isize+1,jsize+1,      0)   =   1
 
      if (( iflags2(      0,      0,ksize  )
     &     +iflags2(      0,      1,ksize+1)
     &     +iflags2(      1,      0,ksize+1)) .eq. 3)
     .      iflags2(      0,      0,ksize+1)   =   1
 
      if (( iflags2(isize+1,      0,ksize  )
     &     +iflags2(isize+1,      1,ksize+1)
     &     +iflags2(isize  ,      0,ksize+1)) .eq. 3)
     .      iflags2(isize+1,      0,ksize+1)   =   1
 
      if (( iflags2(      0,jsize+1,ksize  )
     &     +iflags2(      0,jsize  ,ksize+1)
     &     +iflags2(      1,jsize+1,ksize+1)) .eq. 3)
     .      iflags2(      0,jsize+1,ksize+1)   =   1
 
      if (( iflags2(isize+1,jsize+1,ksize  )
     &     +iflags2(isize+1,jsize  ,ksize+1)
     &     +iflags2(isize  ,jsize+1,ksize+1)) .eq. 3)
     .      iflags2(isize+1,jsize+1,ksize+1)   =   1


      if (dprint) then
         write(outunit,*)" from domup: domflags (after expansion)"
         do 70 kk = 1, ksize
         k = ksize + 1 - kk
         write(outunit,*) 'plane k = ',k
         do 70 jj = 1, jsize
         j = jsize + 1 - jj
         write(outunit,100)(iflags2(i,j,k),i=1,isize)
 100     format(80i1)
 70      continue
      endif

      return
      end
