c
c --------------------------------------------------------------
c
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

      do 10 j = nghost+1, mjtot-nghost
      do 10 i = nghost+1, mitot-nghost
         if (rctold(1,i,j) .ne. goodpt) then
           rectflags(i,j) = badpt
         endif
 10   continue
c
 99   return
      end
