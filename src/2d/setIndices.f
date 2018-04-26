c
c -----------------------------------------------------------------------
c
       subroutine setIndices(ist,iend,jst,jend,ilo,ihi,jlo,jhi,
     &                       ishift,jshift,level)

       use amr_module
       implicit double precision (a-h,o-z)

       dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)

c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!> setIndices:  for periodicity a region that sticks out is wrapped into
!!              another region.
c
c  this is just annoying code this is needed in several places so it became a routine
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


       ist(1) = ilo
       ist(2) = 0
       ist(3) = iregsz(level)

       iend(1) = -1
       iend(2) = iregsz(level)-1
       iend(3) = ihi

       jst(1) = jlo
       jst(2) = 0
       jst(3) = jregsz(level)

       jend(1) = -1
       jend(2) = jregsz(level)-1
       jend(3) = jhi

       if (xperdom) then ! regions that stick out of domain get periodically wrapped indices
          ishift(1) = iregsz(level)
          ishift(2) = 0
          ishift(3) = -iregsz(level)
       else   ! no shifting of indices
          ishift(1) = 0
          ishift(2) = 0
          ishift(3) = 0
       endif

       if (yperdom) then
          jshift(1) = jregsz(level)
          jshift(2) = 0
          jshift(3) = -jregsz(level)
       else  ! no shifting of indices
          jshift(1) = 0
          jshift(2) = 0
          jshift(3) = 0
       endif


       return
       end
