c
c ----------------------------------------------------------
c
       subroutine shiftset(intarray,intarray2,isize,jsize)
c       subroutine old_shiftset(intarray,intarray2,idir,jdir,isize,jsize)

       implicit double precision (a-h, o-z)

       include "call.i"

       integer*1 intarray (0:isize+1,0:jsize+1), 
     1           intarray2(0:isize+1,0:jsize+1)

c :::::::::::::::::::::: CSHIFT :::::::::::::::::::::::::::::::
c shift by + or - 1 in either direction (but only 1 at a time)
c used for bit calculus for proper nesting, buffering, etc.
c similar to cshift on CM machine.
c includes periodic buffering as well.
c
c NEWER VERSION: DO ALL DIRS AT SAME TIME
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       if (xperdom) then
          do 10 j = 0, jsize+1
             intarray(0,j) = intarray(isize,j)
             intarray(isize+1,j) = intarray(1,j)
 10       continue
       else
          do 11 j = 0, jsize+1
             intarray(0,j) = 0
             intarray(isize+1,j) = 0
 11       continue
       endif
       if (yperdom) then
          do 12 i = 0, isize+1
             intarray(i,0) = intarray(i,jsize)
             intarray(i,jsize+1) = intarray(i,1)
 12       continue
       else if (spheredom) then   !use mapped stuff for sphere
          do 14 i = 0, isize+1
             intarray(i,0)       = intarray(isize+1-i,1)
             intarray(i,jsize+1) = intarray(isize+1-i,jsize)
 14       continue
       else
          do 13 i = 0, isize+1
             intarray(i,0) = 0
             intarray(i,jsize+1) = 0
 13       continue
       endif

       do j = 1, jsize
       do i = 1, isize
          intarray2(i,j) = 0
       end do
       end do

       do j = 1, jsize
       do i = 1, isize

          intflag = intarray(i,j)
          if (intflag .gt. 0)  then
c            cell is flagged, buffer in all dirs by one cell
c            use second array to avoid propagation
             mlo = max(i - ibuff, 0)
             mhi = min(i + ibuff, isize+1)
             klo = max(j - ibuff, 0)
             khi = min(j + ibuff, jsize+1)
             do k = klo, khi 
             do m = mlo, mhi
                intarray2(m,k) = intflag  ! copy the flag (may not be = 1?)
             end do
             end do
          endif

       end do
       end do


c       if (idir .eq. 1) then
c           do 22 j = 1, jsize
c           do 22 i = 1, isize
c              intarray2(i,j) = intarray(i+1,j)
c 22        continue
c       elseif (idir .eq. -1) then
c           do 25 j = 1, jsize
c           do 25 i = 1, isize
c               intarray2(i,j) = intarray(i-1,j)
c 25        continue
c       elseif (jdir .eq. 1) then
c           do 50 j = 1, jsize
c           do 50 i = 1, isize
c               intarray2(i,j) = intarray(i,j+1)
c 50         continue
c       elseif (jdir .eq. -1) then
c           do 55 j = 1, jsize
c           do 55 i = 1, isize
c              intarray2(i,j) = intarray(i,j-1)
c 55        continue
c       endif
c

c   copy back. need flags in original array

       do 60 j = 1, jsize
       do 60 i = 1, isize
c        intarray(i,j) = max(intarray(i,j),intarray2(i,j))
         intarray(i,j) = intarray(i,j) + intarray2(i,j)
 60    continue


       return
       end
