c
c ----------------------------------------------------------------
c
       subroutine auxcoarsen(auxdub,midub,mjdub,mkdub,auxbgc,
     1                       mi2tot,mj2tot,mk2tot,naux,auxtype)

       implicit double precision (a-h, o-z)

       dimension     auxdub(naux, midub, mjdub, mkdub)
       dimension     auxbgc(naux, mi2tot,mj2tot,mk2tot)
       character*10  auxtype(naux)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
c coarsen = coarsen the fine grid auxiliary data (with double the usual
c           number of ghost cells to prepare coarsened data
c           for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       do 50 iaux = 1, naux

       if (auxtype(iaux) .eq. "center" .or.
     .     auxtype(iaux) .eq. "capacity") then
         do 20 k = 1, mk2tot
            kfine = 2*(k-1) + 1
            do 20 j = 1, mj2tot
               jfine = 2*(j-1) + 1
               do 20 i = 1, mi2tot
                  ifine = 2*(i-1) + 1
                  auxbgc(iaux,i,j,k) =
     &                      ( auxdub(iaux,ifine  ,jfine  ,kfine  )
     &                       +auxdub(iaux,ifine+1,jfine  ,kfine  )
     &                       +auxdub(iaux,ifine  ,jfine+1,kfine  )
     &                       +auxdub(iaux,ifine+1,jfine+1,kfine  )
     &                       +auxdub(iaux,ifine  ,jfine  ,kfine+1)
     &                       +auxdub(iaux,ifine+1,jfine  ,kfine+1)
     &                       +auxdub(iaux,ifine  ,jfine+1,kfine+1)
     &                       +auxdub(iaux,ifine+1,jfine+1,kfine+1))/8.d0
20       continue

       elseif (auxtype(iaux) .eq. "xleft") then
         do 10 k = 1, mk2tot
            kfine = 2*(k-1) + 1
            do 10 j = 1, mj2tot
               jfine = 2*(j-1) + 1
               do 10 i = 1, mi2tot
                  ifine = 2*(i-1) + 1
                  auxbgc(iaux,i,j,k) =
     &                        ( auxdub(iaux,ifine,jfine  ,kfine  )
     &                         +auxdub(iaux,ifine,jfine+1,kfine  )
     &                         +auxdub(iaux,ifine,jfine  ,kfine+1)
     &                         +auxdub(iaux,ifine,jfine+1,kfine+1))/4.d0
10       continue

       elseif (auxtype(iaux) .eq. "yleft") then
         do 15 k = 1, mk2tot
            kfine = 2*(k-1) + 1
            do 15 j = 1, mj2tot
               jfine = 2*(j-1) + 1
               do 15 i = 1, mi2tot
                  ifine = 2*(i-1) + 1
                  auxbgc(iaux,i,j,k) =
     &                         (auxdub(iaux,ifine  ,jfine,kfine  )
     &                         +auxdub(iaux,ifine+1,jfine,kfine  )
     &                         +auxdub(iaux,ifine  ,jfine,kfine+1)
     &                         +auxdub(iaux,ifine+1,jfine,kfine+1))/4.d0
15       continue

       elseif (auxtype(iaux) .eq. "zleft") then
         do 19 k = 1, mk2tot
            kfine = 2*(k-1) + 1
            do 19 j = 1, mj2tot
               jfine = 2*(j-1) + 1
               do 19 i = 1, mi2tot
                  ifine = 2*(i-1) + 1
                  auxbgc(iaux,i,j,k) =
     &                         (auxdub(iaux,ifine  ,jfine  ,kfine)
     &                         +auxdub(iaux,ifine+1,jfine  ,kfine)
     &                         +auxdub(iaux,ifine  ,jfine+1,kfine)
     &                         +auxdub(iaux,ifine+1,jfine+1,kfine))/4.d0
19       continue

       endif

50     continue

       return
       end
