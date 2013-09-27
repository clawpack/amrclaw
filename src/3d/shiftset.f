c
c ----------------------------------------------------------
c
       subroutine shiftset(intarray,isize,jsize,ksize)

      use amr_module
      implicit double precision (a-h, o-z)

      integer*1 intarray (0:isize+1,0:jsize+1,0:ksize+1)
      integer*1 intarray2(0:isize+1,0:jsize+1,0:ksize+1)

c :::::::::::::::::::::: CSHIFT :::::::::::::::::::::::::::::::
c shift by + or - 1 in either direction (but only 1 at a time)
c used for bit calculus for proper nesting, buffering, etc.
c similar to cshift on CM machine.
c includes periodic buffering as well.
c
c NEWER VERSION: DO ALL DIRS AT SAME TIME
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       if (xperdom) then
          do 10 k = 0, ksize+1
          do 10 j = 0, jsize+1
            intarray(      0,j,k) = intarray(isize,j,k)
            intarray(isize+1,j,k) = intarray(    1,j,k)
 10       continue
       else
          do 11 k = 0, ksize+1
          do 11 j = 0, jsize+1
            intarray(      0,j,k) = 0
            intarray(isize+1,j,k) = 0
 11       continue
       endif
       if (yperdom) then
         do 12 k = 0, ksize+1
         do 12 i = 0, isize+1
            intarray(i,      0,k) = intarray(i,jsize,k)
            intarray(i,jsize+1,k) = intarray(i,    1,k)
 12       continue
       else
         do 13 k = 0, ksize+1
         do 13 i = 0, isize+1
            intarray(i,      0,k) = 0
            intarray(i,jsize+1,k) = 0
 13       continue
       endif
       if (zperdom) then
          do 14 j = 0, jsize+1
          do 14 i = 0, isize+1
             intarray(i,j,      0) = intarray(i,j,ksize)
             intarray(i,j,ksize+1) = intarray(i,j,    1)
 14       continue
       else
          do 15 j = 0, jsize+1
          do 15 i = 0, isize+1
             intarray(i,j,      0) = 0
             intarray(i,j,ksize+1) = 0
 15       continue
       endif

       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize
          intarray2(i,j,k) = 0
       end do
       end do
       end do

       do k = 1, ksize
       do j = 1, jsize
       do i = 1, isize

          intflag = intarray(i,j,k)
          if (intflag .gt. 0)  then
c            cell is flagged, buffer in all dirs by one cell
c            use second array to avoid propagation
             iilo = max(i - ibuff, 0)
             iihi = min(i + ibuff, isize+1)
             jjlo = max(j - ibuff, 0)
             jjhi = min(j + ibuff, jsize+1)
             kklo = max(k - ibuff, 0)
             kkhi = min(k + ibuff, ksize+1)
             do kk = kklo, kkhi 
             do jj = jjlo, jjhi
             do ii = iilo, iihi
                intarray2(ii,jj,kk) = intflag  ! copy the flag (may not be = 1?)
             end do
             end do
             end do
          endif

       end do
       end do
       end do

c   copy back. need flags in original array

       do 60 k = 1, ksize
       do 60 j = 1, jsize
       do 60 i = 1, isize
c         intarray(i,j,k) = max(intarray(i,j,k),intarray2(i,j,k))
         intarray(i,j,k) = intarray(i,j,k) + intarray2(i,j,k)
 60    continue


       return
       end
