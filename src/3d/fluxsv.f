c
c ----------------------------------------------------------
c
      subroutine fluxsv(mptr,xfluxm,xfluxp,yfluxm,yfluxp,zfluxm,zfluxp,
     1                  listbc,
     2                  ndimx,ndimy,ndimz,nvar,maxsp,dtc,hx,hy,hz)
c
      use amr_module
      implicit double precision (a-h,o-z)

      parameter(numbcs=6)

      dimension xfluxm(nvar,ndimx,ndimy,ndimz)
      dimension xfluxp(nvar,ndimx,ndimy,ndimz)
      dimension yfluxm(nvar,ndimx,ndimy,ndimz)
      dimension yfluxp(nvar,ndimx,ndimy,ndimz)
      dimension zfluxm(nvar,ndimx,ndimy,ndimz)
      dimension zfluxp(nvar,ndimx,ndimy,ndimz)
      dimension listbc(numbcs,maxsp)
c
c :::::::::::::::::::: FLUXSV :::::::::::::::::::::::::
c
c  coarse grids should save their fluxes in cells adjacent to
c  their nested fine grids, for later conservation fixing.
c  listbc holds info for where to save which fluxes.
c  xflux holds 'f' fluxes, yflux holds 'g' fluxes, zflux holds 'h' fluxes.
c
c :::::::::::::::::::::::::::::;:::::::::::::::::::::::
 
 
      ispot   = 1
      level   = node(nestlevel,mptr)

 10      if (listbc(1,ispot).eq.0) go to 99          
c
         mkid     = listbc(5,ispot)
         intopl   = listbc(6,ispot)
         nx       = node(ndihi,mkid) - node(ndilo,mkid) + 1
         ny       = node(ndjhi,mkid) - node(ndjlo,mkid) + 1
         nz       = node(ndkhi,mkid) - node(ndklo,mkid) + 1
         kidlst   = node(ffluxptr,mkid)
         i        = listbc(1,ispot)
         j        = listbc(2,ispot)
         k        = listbc(3,ispot)
         inlist   = kidlst + nvar*(intopl-1) - 1

         if (listbc(4,ispot) .eq. 1) then
c           ::::: Cell i,j,k is on right side of a fine grid
            do 100 ivar = 1, nvar
              alloc(inlist + ivar) = -xfluxp(ivar,i  ,j  ,k  )*dtc*hy*hz
100         continue
         endif

         if (listbc(4,ispot) .eq. 2) then
c           ::::: Cell i,j,k is on front side of fine grid
            do 200 ivar = 1, nvar
              alloc(inlist + ivar) = -yfluxm(ivar,i  ,j+1,k  )*dtc*hx*hz
200         continue
         endif

         if (listbc(4,ispot) .eq. 3) then
c           ::::: Cell i,j,k is on left side of fine grid
            do 300 ivar = 1, nvar
              alloc(inlist + ivar) = -xfluxm(ivar,i+1,j  ,k  )*dtc*hy*hz
300         continue
         endif

         if (listbc(4,ispot) .eq. 4) then
c           ::::: Cell i,j,k is on rear side of fine grid
            do 400 ivar = 1, nvar
              alloc(inlist + ivar) = -yfluxp(ivar,i  ,j  ,k  )*dtc*hx*hz
400         continue
         endif

         if (listbc(4,ispot) .eq. 5) then
c           ::::: Cell i,j,k is on top side of fine grid
            do 500 ivar = 1, nvar
              alloc(inlist + ivar) = -zfluxp(ivar,i  ,j  ,k  )*dtc*hx*hy
500         continue
         endif

         if (listbc(4,ispot) .eq. 6) then
c           ::::: Cell i,j,k is on bottom side of fine grid
            do 600 ivar = 1, nvar
              alloc(inlist + ivar) = -zfluxm(ivar,i  ,j  ,k+1)*dtc*hx*hy
600         continue
         endif
c
      ispot = ispot + 1
      if (ispot .gt. maxsp) go to 99
      go to 10
c
 99   return
      end
