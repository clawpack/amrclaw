c
c ----------------------------------------------------------
c
      subroutine fluxsv(mptr,xfluxm,xfluxp,listbc,
     1                  ndimx,nvar,maxsp,dtc,hx)
c
      use amr_module
      implicit double precision (a-h,o-z)


      dimension xfluxp(nvar,ndimx)
      dimension xfluxm(nvar,ndimx)
      dimension listbc(4,maxsp)
c
c :::::::::::::::::::: FLUXSV :::::::::::::::::::::::::
c
c  coarse grids should save their fluxes in cells adjacent to
c  their nested fine grids, for later conservation fixing.
c  listbc holds info for where to save which fluxes.
c  xflux holds 'f' fluxes.
c
c :::::::::::::::::::::::::::::;:::::::::::::::::::::::
 
 
      ispot   = 1
      level   = node(nestlevel,mptr)

 10      if (listbc(1,ispot).eq.0) go to 99          
c
         mkid     = listbc(3,ispot)
         intopl   = listbc(5,ispot)
         nx       = node(ndihi,mkid) - node(ndilo,mkid) + 1
         kidlst   = node(ffluxptr,mkid)
         i        = listbc(1,ispot)
         inlist   = kidlst + nvar*(intopl-1) - 1
c
c side j (listbc 2) has which side of coarse cell has interface
c so can save appropriate fluxes.  (dont know why we didnt have
c which flux to save directly (i.e. put i+1,j to save that flux
c rather than putting in cell center coords).

         if (listbc(2,ispot) .eq. 1) then
c           ::::: Cell i is on right side of a fine grid
            do 100 ivar = 1, nvar
               alloc(inlist + ivar) = -xfluxp(ivar,i)*dtc
100         continue
c         write(dbugunit,901) i,j,1,(xfluxp(ivar,i,j),ivar=1,nvar)
         endif

         if (listbc(2,ispot) .eq. 2) then
c           ::::: Cell i on left side of fine grid
            do 300 ivar = 1, nvar
               alloc(inlist + ivar) = -xfluxm(ivar,i+1)*dtc
300         continue
c         write(dbugunit,901) i,j,3,(xfluxm(ivar,i+1,j),ivar=1,nvar)
         endif
c
c        ### new bcs 5 and 6 come from spherical mapping. note sign change:
c        ### previous fluxes stored negative flux, fine grids always add
c        ### their flux, then the delta is either added or subtracted as
c        ### appropriate for that side.  New bc adds or subtracts BOTH fluxes.
c
         if (listbc(2,ispot) .eq. 5) then
c           ::::: Cell i,j on top side of fine grid with spherical mapped bc
c            do 500 ivar = 1, nvar
c               alloc(inlist + ivar) = yfluxm(ivar,i,j+1)*dtc*hx
c500         continue
c         write(dbugunit,901) i,j,5,(yfluxm(ivar,i,j+1),ivar=1,nvar)
c 901        format(2i4," side",i3,4e15.7)
             write(*,*) "Error, spherical mapping not enabled yet."
         endif
c
         if (listbc(2,ispot) .eq. 6) then
c           ::::: Cell i,j on bottom side of fine grid with spherical mapped bc
c            do 600 ivar = 1, nvar
c               alloc(inlist + ivar) = yfluxp(ivar,i,j)*dtc*hx
c600         continue
c         write(dbugunit,901) i,j,6,(yfluxp(ivar,i,j),ivar=1,nvar)
             write(*,*) "Error, spherical mapping not enabled yet."
         endif

      ispot = ispot + 1
      if (ispot .gt. maxsp) go to 99
      go to 10
c
 99   return
      end
