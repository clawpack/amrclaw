c
!> Transfer flagged arrays from errest into the one from spatial differencing 
!!
!! NOTE: not dimensioned the same. rectflags is possibly larger to accomodate
!! in-place buffering.

      subroutine addflags(rectflags,mibuff,mjbuff,
     1                    rctold,idim3,mitot,mjtot,mptr)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension rctold(idim3,mitot,mjtot)
      dimension rectflags(mibuff,mjbuff)

c :::::::::::::::::::::: ADDFLAGS ::::::::::::::::::::::::::::::::::
c transfer flagged arrays from errest into the one from spatial
c differencing 
c NOTE: not dimensioned the same. rectflags is possibly larger to accomodate
c in-place buffering.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      do j = nghost+1, mjtot-nghost
      do i = nghost+1, mitot-nghost
         if (rctold(1,i,j) .gt. DONTFLAG) then
           rectflags(i,j) = DOFLAG
         endif
      end do
      end do
c
 99   return
      end
