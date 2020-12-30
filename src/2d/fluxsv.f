!> When a coarse grid cell is advanced, if it borders a coarse-fine
!! interface, the flux or wave that emits from the interface and goes into the coarse cell 
!! is added to the corresponding location in **node(ffluxptr, mkid)**
!! for conservative fix later, where **mkid** is grid number of the
!! nested fine grid on the other side of that interface. 
!! 
!! 
c ----------------------------------------------------------
c
      subroutine fluxsv(mptr,xfluxm,xfluxp,yfluxm,yfluxp,listbc,
     1                  ndimx,ndimy,nvar,maxsp,dtc,hx,hy)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension xfluxp(nvar,ndimx,ndimy), yfluxp(nvar,ndimx,ndimy)
      dimension xfluxm(nvar,ndimx,ndimy), yfluxm(nvar,ndimx,ndimy)
      !dimension listbc(5,maxsp)
      double precision listbc(5,maxsp)
c
c :::::::::::::::::::: FLUXSV :::::::::::::::::::::::::
c
c  coarse grids should save their fluxes in cells adjacent to
c  their nested fine grids, for later conservation fixing.
c  listbc holds info for where to save which fluxes.
c  xflux holds 'f' fluxes, yflux holds 'g' fluxes.
c
c :::::::::::::::::::::::::::::;:::::::::::::::::::::::
 
 
      ispot   = 1
      level   = node(nestlevel,mptr)

 10      if (listbc(1,ispot).eq.0) go to 99          
c
         mkid     = listbc(4,ispot)
         intopl   = listbc(5,ispot)
         nx       = node(ndihi,mkid) - node(ndilo,mkid) + 1
         ny       = node(ndjhi,mkid) - node(ndjlo,mkid) + 1
         kidlst   = node(ffluxptr,mkid)
         i        = listbc(1,ispot)
         j        = listbc(2,ispot)
         inlist   = kidlst + nvar*(intopl-1) - 1
c
c side k (listbc 3) has which side of coarse cell has interface
c so can save appropriate fluxes.  (dont know why we didnt have
c which flux to save directly (i.e. put i+1,j to save that flux
c rather than putting in cell center coords).

         if (listbc(3,ispot) .eq. 1) then
c           ::::: Cell i,j is on right side of a fine grid
            do 100 ivar = 1, nvar
               alloc(inlist + ivar) = -xfluxp(ivar,i,j)*dtc*hy
100         continue
c         write(dbugunit,901) i,j,1,(xfluxp(ivar,i,j),ivar=1,nvar)
         endif

         if (listbc(3,ispot) .eq. 2) then
c           ::::: Cell i,j on bottom side of fine grid
            do 200 ivar = 1, nvar
               alloc(inlist + ivar) = -yfluxm(ivar,i,j+1)*dtc*hx
200         continue
c         write(dbugunit,901) i,j,2,(yfluxm(ivar,i,j+1),ivar=1,nvar)
         endif

         if (listbc(3,ispot) .eq. 3) then
c           ::::: Cell i,j on left side of fine grid
            do 300 ivar = 1, nvar
               alloc(inlist + ivar) = -xfluxm(ivar,i+1,j)*dtc*hy
300         continue
c         write(dbugunit,901) i,j,3,(xfluxm(ivar,i+1,j),ivar=1,nvar)
         endif

         if (listbc(3,ispot) .eq. 4) then
c           ::::: Cell i,j on top side of fine grid
            do 400 ivar = 1, nvar
               alloc(inlist + ivar) = -yfluxp(ivar,i,j)*dtc*hx
400         continue
c         write(dbugunit,901) i,j,4,(yfluxp(ivar,i,j),ivar=1,nvar)
         endif
c
c        ### new bcs 5 and 6 come from spherical mapping. note sign change:
c        ### previous fluxes stored negative flux, fine grids always add
c        ### their flux, then the delta is either added or subtracted as
c        ### appropriate for that side.  New bc adds or subtracts BOTH fluxes.
c
         if (listbc(3,ispot) .eq. 5) then
c           ::::: Cell i,j on top side of fine grid with spherical mapped bc
            do 500 ivar = 1, nvar
               alloc(inlist + ivar) = yfluxm(ivar,i,j+1)*dtc*hx
500         continue
c         write(dbugunit,901) i,j,5,(yfluxm(ivar,i,j+1),ivar=1,nvar)
 901        format(2i4," side",i3,4e15.7)
         endif
c
         if (listbc(3,ispot) .eq. 6) then
c           ::::: Cell i,j on bottom side of fine grid with spherical mapped bc
            do 600 ivar = 1, nvar
               alloc(inlist + ivar) = yfluxp(ivar,i,j)*dtc*hx
600         continue
c         write(dbugunit,901) i,j,6,(yfluxp(ivar,i,j),ivar=1,nvar)
         endif

      ispot = ispot + 1
      if (ispot .gt. maxsp) go to 99
      go to 10
c
 99   return
      end
