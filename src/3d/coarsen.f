


c
c ---------------------------------------------------
c
       subroutine coarsen(valdub,midub,mjdub,mkdub,
     &                    valbgc,mi2tot,mj2tot,mk2tot,nvar)
      
       implicit double precision (a-h, o-z)

       dimension valdub(midub ,mjdub ,mkdub ,nvar)
       dimension valbgc(mi2tot,mj2tot,mk2tot,nvar)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
c coarsen = coarsen the fine grid data (with double the usual
c           number of ghost cells to prepare coarsened
c           grid for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       do 10 ivar = 1, nvar
       do 10 k = 1, mk2tot
       kfine = 2*(k-1) + 1
       do 10 j = 1, mj2tot
       jfine = 2*(j-1) + 1
       do 10 i = 1, mi2tot
       ifine = 2*(i-1) + 1
       valbgc(i,j,k,ivar) = (valdub(ifine  ,jfine  ,kfine  ,ivar)
     &                      +valdub(ifine+1,jfine  ,kfine  ,ivar)
     &                      +valdub(ifine  ,jfine+1,kfine  ,ivar)
     &                      +valdub(ifine+1,jfine+1,kfine  ,ivar)
     &                      +valdub(ifine  ,jfine  ,kfine+1,ivar)
     &                      +valdub(ifine+1,jfine  ,kfine+1,ivar)
     &                      +valdub(ifine  ,jfine+1,kfine+1,ivar)
     &                      +valdub(ifine+1,jfine+1,kfine+1,ivar))/8.d0
10     continue

       return
       end
