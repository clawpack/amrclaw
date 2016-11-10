c
c ----------------------------------------------------------------
c
       subroutine auxcoarsen(auxdub,midub,auxbgc,
     1                       mi2tot,naux,auxtype)
      
       implicit double precision (a-h, o-z)

       dimension     auxdub(naux,midub)
       dimension     auxbgc(naux,mi2tot)
       character*10  auxtype(naux)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
c coarsen = coarsen the fine grid auxiliary data (with double the usual
c           number of ghost cells to prepare coarsened data
c           for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       do 50 iaux = 1, naux

       if (auxtype(iaux) .eq. "center" .or. 
     .     auxtype(iaux) .eq. "capacity") then
            do 20 i = 1, mi2tot
               ifine = 2*(i-1) + 1
               auxbgc(iaux,i) = (auxdub(iaux,ifine) +
     &                             auxdub(iaux,ifine+1))/2.d0
20       continue

       elseif (auxtype(iaux) .eq. "xleft") then
            do 10 i = 1, mi2tot
               ifine = 2*(i-1) + 1
               auxbgc(iaux,i) = auxdub(iaux,ifine)
10       continue

       endif

50     continue

       return
       end
