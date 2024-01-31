c
c ----------------------------------------------------
c
      subroutine domup(iflags2,iflags,ibase,jbase,isize,jsize,lev)

      use amr_module
      implicit double precision (a-h, o-z)


      integer(kind=1)  iflags2(0:isize+1,0:jsize+1)
      integer(kind=1)  iflags (0:ibase+1,0:jbase+1)

c
c ::::::::::::::::::::::::::: DOMUP :::::::::::::::::::::
c 
!>  domain flags are in iflags. copy into iflags2, allowing
!!  for change of level and dimension
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::

      if (dprint) then
         write(outunit,*)" from domup: domflags (before expansion)"
         do 5 jj = 1, jbase
         j = jbase + 1 - jj
         write(outunit,100)(iflags(i,j),i=1,ibase)
 5       continue
      endif

      do j = 0, jsize+1
      do i = 0, isize+1
         iflags2(i,j) = 0
      end do
      end do

      do j = 1, jbase
      do i = 1, ibase
          ifine = (i-1) * intratx(lev)
          jfine = (j-1) * intraty(lev)
          do mj = 1, intraty(lev)
          do mi = 1, intratx(lev)
            iflags2(ifine+mi,jfine+mj) = iflags(i,j)  
          end do
          end do
      end do
      end do
c
c  take care of periodicity again or if border of domain touches a 
c  physical boundary then set domain in ghost cell as well
c
      if (xperdom) then
         do 35 j = 0, jsize+1
           iflags2(0,j)       = iflags2(isize,j)
           iflags2(isize+1,j) = iflags2(1,j)
 35      continue
       else
       do 55 j = 1, jsize
         if (iflags2(1,j) .eq. 1) iflags2(0,j) = 1
         if (iflags2(isize,j) .eq. 1) iflags2(isize+1,j) = 1
 55    continue
      endif
      if (yperdom) then
         do 45 i = 0, isize+1
           iflags2(i,0)       = iflags2(i,jsize)
           iflags2(i,jsize+1) = iflags2(i,1)
 45      continue
       else if (spheredom) then
         do 46 i = 0, isize+1
           iflags2(i,0)       = iflags2(isize+1-i,1)
           iflags2(i,jsize+1) = iflags2(isize+1-i,jsize)
 46      continue

       else
         do 65 i = 1, isize
           if (iflags2(i,1) .eq. 1) iflags2(i,0) = 1
           if (iflags2(i,jsize) .eq. 1) iflags2(i,jsize+1) = 1
 65      continue
      endif

c
c the 4 corners
c
        if (iflags2(0,1)+iflags2(1,0) .eq. 2) iflags2(0,0)=1
        if (iflags2(isize,0)+iflags2(isize+1,1) .eq. 2)
     .          iflags2(isize+1,0)=1
        if (iflags2(isize,jsize+1)+iflags2(isize+1,jsize) .eq. 2)
     .          iflags2(isize+1,jsize+1)=1
        if (iflags2(0,jsize)+iflags2(1,jsize+1) .eq. 2)
     .          iflags2(0,jsize+1)=1


      if (dprint) then
         write(outunit,*)" from domup: domflags (after expansion)"
         do 70 jj = 1, jsize
         j = jsize + 1 - jj
         write(outunit,100)(iflags2(i,j),i=1,isize)
 100     format(80i1)
 70      continue
      endif

      return
      end
