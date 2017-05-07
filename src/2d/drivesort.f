!> Sort 2D points (stored in badpts) based on their equivalent 
!! integer value (based on their x,y coordinates). Also remove
!! duplicate points in badpts.
c
c -------------------------------------------------------------
c
       subroutine drivesort(npts,badpts,level,index,mbuff)

      use amr_module
      implicit  double precision (a-h,o-z)
      dimension badpts(2,npts)
      dimension iflags(npts), ixArray(npts)
      logical db/.false./
 
      iadd(i,j) = (i+mbuff)  + (isize+2*mbuff)*(j+mbuff)
c
c  convert using one dimensional ordering of badpts array as if
c  it covered entire domain (xprob by yprob) on this level
c
      isize = iregsz(level)
      jsize = jregsz(level)

      do k = 1, npts
        i = badpts(1,k)-.5  ! remember was shifted when put  into array
        j = badpts(2,k)-.5  
        intEquiv = iadd(i,j)
c        write(*,*)i,j," has equivalent integer ",intEquiv
        iflags(k) = intEquiv
      end do

      call qsorti(ixArray, npts, iflags)

c copy back to badpts, in sorted order, removing duplicates      
      k = 1
      index = 0
      do while (k .le. npts) 
         intEquiv = iflags(ixArray(k))
         index = index + 1
         badpts(2,index) = intEquiv/(isize+2*mbuff) + .5 -mbuff
         badpts(1,index) = mod(intEquiv,(isize+2*mbuff)) + .5 -mbuff
        if (db) write(outunit,101) badpts(1,index),badpts(2,index)
 101    format(2f6.1)
         k = k + 1
         do while ( k.le. npts)    ! skip over duplicates
            if (iflags(ixArray(k)) .eq. iflags(ixArray(k-1))) then
c           write(*,*)" duplicate in sorted array loc ",ixarray(k)
               k = k+1
            else
               exit   ! back to outer loop
            endif
         end do
       if (k .gt. npts) exit !did we stop because we ran off end or pts not equal
      end do

      if (gprint) then
          write(outunit,929) index
 929      format(i5," flagged pts after removing duplicates and ",
     &           " non-nested flags")
      endif

      return
      end
