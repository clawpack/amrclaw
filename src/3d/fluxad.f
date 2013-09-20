c
c -------------------------------------------------------
c
      subroutine fluxad(xfluxm,xfluxp,yfluxm,yfluxp,zfluxm,zfluxp,
     1                  svdflx,mptr,mitot,mjtot,mktot,
     2                  nvar,lenbc,lratiox,lratioy,lratioz,
     3                  ng,dtf,dx,dy,dz)
c

      use amr_module
      implicit double precision (a-h,o-z)

c :::::::::::::::::::: FLUXAD ::::::::::::::::::::::::::::::::::
c  save fine grid fluxes  at the border of the grid, for fixing
c  up the adjacent coarse cells. at each edge of the grid, only
c  save the plus or minus fluxes, as necessary. For ex., on
c  left edge of fine grid, it is the minus xfluxes that modify the
c  coarse cell.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      dimension xfluxm(nvar,mitot,mjtot,mktot)
      dimension xfluxp(nvar,mitot,mjtot,mktot)
      dimension yfluxm(nvar,mitot,mjtot,mktot)
      dimension yfluxp(nvar,mitot,mjtot,mktot)
      dimension zfluxm(nvar,mitot,mjtot,mktot)
      dimension zfluxp(nvar,mitot,mjtot,mktot)
      dimension svdflx(nvar,lenbc)
 
      nx  = mitot-2*ng
      ny  = mjtot-2*ng
      nz  = mktot-2*ng
      nzc = nz/lratioz
      nyc = ny/lratioy
      nxc = nx/lratiox
 
c ::::: left side saved first
      lind = 0

      do 100 k=1,nzc
         kfine = (k-1)*lratioz + ng
         do 110 j=1,nyc
            jfine = (j-1)*lratioy + ng
            lind = lind + 1
            do 120 ivar = 1, nvar
               do 130 kl=1,lratioz
               do 130 jl=1,lratioy
                  svdflx(ivar,lind) = svdflx(ivar,lind) +
     1                  xfluxm(ivar,ng+1,jfine+jl,kfine+kl)*dtf*dy*dz
c	          write(dbugunit,900)lind,
c     .                              xfluxm(ivar,1,jfine+jl,kfine+kl),
c     .                              xfluxp(ivar,1,jfine+jl,kfine+kl)
130            continue
120         continue
110      continue
100   continue
 
c ::::: rear side
      do 200 k=1,nzc
         kfine = (k-1)*lratioz + ng
         do 210 i=1,nxc
            ifine = (i-1)*lratiox + ng
            lind = lind + 1
            do 220 ivar = 1, nvar
               do 230 kl=1,lratioz
               do 230 il=1,lratiox
                  svdflx(ivar,lind) = svdflx(ivar,lind) + 
     1              yfluxp(ivar,ifine+il,mjtot-ng+1,kfine+kl)*dtf*dx*dz
c         	   write(dbugunit,900)lind,
c     .                       yfluxm(ivar,ifine+il,mjtot-ng+1,kfine+kl),
c     .                       yfluxp(ivar,ifine+il,mjtot-ng+1,kfine+kl)
230            continue
220         continue
210      continue
200   continue
 
c ::::: right side
      do 300 k=1,nzc
         kfine = (k-1)*lratioz + ng
         do 310 j=1,nyc
            jfine = (j-1)*lratioy + ng
            lind = lind + 1
            do 320 ivar = 1, nvar
               do 330 kl=1,lratioz
               do 330 jl=1,lratioy
                  svdflx(ivar,lind) = svdflx(ivar,lind) + 
     1              xfluxp(ivar,mitot-ng+1,jfine+jl,kfine+kl)*dtf*dy*dz
c                  write(dbugunit,900)lind
c     .                       xfluxm(ivar,mitot-ng+1,jfine+jl,kfine+kl),
c     .                       xfluxp(ivar,mitot-ng+1,jfine+jl,kfine+kl)
330            continue
320         continue
310      continue
300   continue
 
c ::::: front side
      do 400 k=1,nzc
         kfine = (k-1)*lratioz + ng
         do 410 i=1,nxc
            ifine = (i-1)*lratiox + ng
            lind = lind + 1
            do 420 ivar = 1, nvar
               do 430 kl=1,lratioz
               do 430 il=1,lratiox
                  svdflx(ivar,lind) = svdflx(ivar,lind) +
     1                    yfluxm(ivar,ifine+il,ng+1,kfine+kl)*dtf*dx*dz
c                  write(dbugunit,900)lind,
c     .                             yfluxm(ivar,ifine+il,ng+1,kfine+kl),
c     .                             yfluxp(ivar,ifine+il,ng+1,kfine+kl)
430            continue
420         continue
410      continue
400   continue
 
c ::::: bottom side
      do 500 j=1,nyc
         jfine = (j-1)*lratioy + ng
         do 510 i=1,nxc
            ifine = (i-1)*lratiox + ng
            lind = lind + 1
            do 520 ivar = 1, nvar
               do 530 jl=1,lratioy
               do 530 il=1,lratiox
                  svdflx(ivar,lind) = svdflx(ivar,lind) +
     1                    zfluxm(ivar,ifine+il,jfine+jl,ng+1)*dtf*dx*dy
c                  write(dbugunit,900)lind,
c     .                             zfluxm(ivar,ifine+il,jfine+jl,ng+1),
c     .                             zfluxp(ivar,ifine+il,jfine+jl,ng+1) 
530            continue
520         continue
510      continue
500   continue
 
c ::::: top side
      do 600 j=1,nyc
         jfine = (j-1)*lratioy + ng
         do 610 i=1,nxc
            ifine = (i-1)*lratiox + ng
            lind = lind + 1
            do 620 ivar = 1, nvar
               do 630 jl=1,lratioy
               do 630 il=1,lratiox
                  svdflx(ivar,lind) = svdflx(ivar,lind) +
     1              zfluxp(ivar,ifine+il,jfine+jl,mktot-ng+1)*dtf*dx*dy
c                  write(dbugunit,900)lind,
c     .                             zfluxm(ivar,ifine+il,jfine+jl,ng+1),
c     .                             zfluxp(ivar,ifine+il,jfine+jl,ng+1)
630            continue
620         continue
610      continue
600   continue

 900              format(' lind ', i4,' m & p ',2e15.7)

      return
      end
