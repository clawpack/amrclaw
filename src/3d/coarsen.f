c
c ---------------------------------------------------
c
       subroutine coarsen(valdub,midub,mjdub,mkdub,
     &                    valbgc,mi2tot,mj2tot,mk2tot,nvar)
      
       implicit double precision (a-h, o-z)

       dimension valdub(nvar,midub ,mjdub ,mkdub )
       dimension valbgc(nvar,mi2tot,mj2tot,mk2tot)

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
       valbgc(ivar,i,j,k) = (valdub(ivar,ifine  ,jfine  ,kfine  )
     &                      +valdub(ivar,ifine+1,jfine  ,kfine  )
     &                      +valdub(ivar,ifine  ,jfine+1,kfine  )
     &                      +valdub(ivar,ifine+1,jfine+1,kfine  )
     &                      +valdub(ivar,ifine  ,jfine  ,kfine+1)
     &                      +valdub(ivar,ifine+1,jfine  ,kfine+1)
     &                      +valdub(ivar,ifine  ,jfine+1,kfine+1)
     &                      +valdub(ivar,ifine+1,jfine+1,kfine+1))/8.d0
10     continue

       return
       end
