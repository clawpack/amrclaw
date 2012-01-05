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
c NOTE: not dimensioned the same. rectflags is larger to accomodate
c inplace buffering
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      mbuff = max(nghost, ibuff) - nghost
      do 10 j = nghost+1, mjtot-nghost
      do 10 i = nghost+1, mitot-nghost
        ioff = mbuff
        rectflags(i+ioff,j+joff) = rectflags(i,j) + rctold(1,i,j)
 10   continue
c
 99   return
      end
