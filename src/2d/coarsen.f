c
c ---------------------------------------------------
c
       subroutine coarsen(valdub,midub,mjdub,valbgc,mi2tot,mj2tot,nvar)
      
       implicit double precision (a-h, o-z)

       dimension  valdub(nvar,midub, mjdub)
       dimension  valbgc(nvar,mi2tot,mj2tot)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
!> coarsen the fine grid data (with double the usual
!! number of ghost cells to prepare coarsened
!! grid for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       do 10 j = 1, mj2tot

          jfine = 2*(j-1) + 1

          do 10 i = 1, mi2tot
            ifine = 2*(i-1) + 1

            do 10 ivar = 1, nvar
  
             valbgc(ivar,i,j) = (valdub(ivar,ifine,jfine) +
     &                           valdub(ivar,ifine+1,jfine)+
     &                           valdub(ivar,ifine,jfine+1) +
     &                           valdub(ivar,ifine+1,jfine+1))/4.d0
10     continue

       return
       end
