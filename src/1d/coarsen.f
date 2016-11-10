c
c ---------------------------------------------------
c
       subroutine coarsen(valdub,midub,valbgc,mi2tot,nvar)
      
       implicit double precision (a-h, o-z)

       dimension  valdub(nvar,midub)
       dimension  valbgc(nvar,mi2tot)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
c coarsen = coarsen the fine grid data (with double the usual
c           number of ghost cells to prepare coarsened
c           grid for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
          do 10 i = 1, mi2tot
            ifine = 2*(i-1) + 1

            do 10 ivar = 1, nvar
  
             valbgc(ivar,i) = (valdub(ivar,ifine) +
     &                           valdub(ivar,ifine+1))/2.d0
10     continue

       return
       end
