c
c ----------------------------------------------------------------
c
       subroutine auxcoarsen(auxdub,midub,mjdub,auxbgc,
     1                       mi2tot,mj2tot,naux,auxtype)
      
       implicit double precision (a-h, o-z)

       dimension     auxdub(naux,midub, mjdub)
       dimension     auxbgc(naux,mi2tot,mj2tot)
       character*10  auxtype(naux)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
!> Coarsen the fine grid auxiliary data (with double the usual
!! number of ghost cells to prepare coarsened data
!! for error estimation).
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       do iaux = 1, naux

       if (auxtype(iaux) .eq. "center" .or. 
     .     auxtype(iaux) .eq. "capacity") then
         do j = 1, mj2tot
            jfine = 2*(j-1) + 1
            do i = 1, mi2tot
               ifine = 2*(i-1) + 1
               auxbgc(iaux,i,j) = (auxdub(iaux,ifine,jfine) +
     &                             auxdub(iaux,ifine+1,jfine)+
     &                             auxdub(iaux,ifine,jfine+1) +
     &                             auxdub(iaux,ifine+1,jfine+1))/4.d0
            end do
         end do

       elseif (auxtype(iaux) .eq. "xleft") then 
         do j = 1, mj2tot
            jfine = 2*(j-1) + 1
            do i = 1, mi2tot
               ifine = 2*(i-1) + 1
               auxbgc(iaux,i,j) = (auxdub(iaux,ifine,jfine) +
     &                             auxdub(iaux,ifine,jfine+1)) /2.d0 
            end do
         end do

       elseif (auxtype(iaux) .eq. "yleft") then 
         do j = 1, mj2tot
            jfine = 2*(j-1) + 1
            do i = 1, mi2tot
               ifine = 2*(i-1) + 1
               auxbgc(iaux,i,j) = (auxdub(iaux,ifine,jfine) +
     &                             auxdub(iaux,ifine+1,jfine))/2.d0
            end do
         end do

       endif

       end do

       return
       end
