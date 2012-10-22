c
c ----------------------------------------------------------
c
       subroutine shiftset(intarray,intarray2,idir ,jdir ,kdir ,
     &                                        isize,jsize,ksize)

       implicit double precision (a-h, o-z)

       include "call.i"

       dimension intarray (0:isize+1,0:jsize+1,0:ksize+1), 
     1		 intarray2(0:isize+1,0:jsize+1,0:ksize+1)

c :::::::::::::::::::::: CSHIFT :::::::::::::::::::::::::::::::
c shift by + or - 1 in either direction (but only 1 at a time)
c used for bit calculus for proper nesting, buffering, etc.
c similar to cshift on CM machine.
c includes periodic buffering as well.
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

       if (idir .eq. 1) then
           do 22 k = 1, ksize
	   do 22 j = 1, jsize
	   do 22 i = 1, isize
	      intarray2(i,j,k) = intarray(i+1,j,k)
 22        continue
       elseif (idir .eq. -1) then
           do 25 k = 1, ksize
	   do 25 j = 1, jsize
	   do 25 i = 1, isize
	      intarray2(i,j,k) = intarray(i-1,j,k)
 25        continue
       elseif (jdir .eq. 1) then
           do 50 k = 1, ksize
	   do 50 j = 1, jsize
	   do 50 i = 1, isize
	      intarray2(i,j,k) = intarray(i,j+1,k)
 50         continue
       elseif (jdir .eq. -1) then
           do 55 k = 1, ksize
	   do 55 j = 1, jsize
	   do 55 i = 1, isize
	      intarray2(i,j,k) = intarray(i,j-1,k)
 55        continue
       elseif (kdir .eq.  1) then
           do 57 k = 1, ksize
           do 57 j = 1, jsize
           do 57 i = 1, isize
              intarray2(i,j,k) = intarray(i,j,k+1)
 57         continue
       elseif (kdir .eq. -1) then
           do 59 k = 1, ksize
           do 59 j = 1, jsize
           do 59 i = 1, isize
              intarray2(i,j,k) = intarray(i,j,k-1)
 59         continue
       endif

c   copy back

       do 60 k = 1, ksize
       do 60 j = 1, jsize
       do 60 i = 1, isize
         intarray(i,j,k) = max(intarray(i,j,k),intarray2(i,j,k))
 60    continue


       return
       end
