c
c ------------------------------------------------------------
c
       subroutine upbnd(listbc,val,nvar,maux,mitot,mjtot,mktot,
     1                  maxsp,mptr)

      use amr_module
      implicit double precision (a-h,o-z)

      parameter(numbcs=6)

       dimension val(nvar,mitot,mjtot,mktot),listbc(numbcs,maxsp),
     1           iused(mitot,mjtot,mktot)
       dimension chsign(numbcs)
       data      chsign/ 1.,-1.,-1., 1., 1.,-1./

       iaddaux(i,j,k) = locaux +     (mcapa-1)
     &                         +     (i-1)*maux
     &                         +     (j-1)*maux*mitot
     &                         +     (k-1)*maux*mitot*mjtot

c
c :::::::::::::::::::::::::::: UPBND :::::::::::::::::::::::::::::
c We now correct the coarse grid with the flux differences stored
c with each of the fine grids. We use an array   iused
c to indicate whether the flux has been updated or not for that zone.
c iused(i,j) = sum from (l=1,4) i(l)*2**(l-1), where i(l) = 1 if the
c flux for the  l-th side of the (i,j)-th cell has already been
c updated, and i(l) = 0 if not.

c if there is a capacity fn. it needs to be included in update formula
c indicated by mcapa not zero (is index of capacity fn.)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      do 10 k=1,mktot
      do 10 j=1,mjtot
      do 10 i=1,mitot
         iused(i,j,k) = 0.
 10   continue

      locaux = node(storeaux,mptr)
      levc   = node(nestlevel,mptr)
      volume = hxposs(levc)*hyposs(levc)*hzposs(levc)
      if (uprint) write(outunit,*)" upbounding grid ",mptr

      do 40 ispot = 1,maxsp
         icrse = listbc(1,ispot)
         if (icrse.eq.0) go to 99

         jcrse = listbc(2,ispot)
         kcrse = listbc(3,ispot)
         iside = listbc(4,ispot)
         norm = 2**(iside-1)
         iflag =iused(icrse,jcrse,kcrse)/norm
         if (mod(iflag,2).eq.1) then
           go to 40
         endif
         mkid = listbc(5,ispot)
         lkid = listbc(6,ispot)
         sgnm = chsign(iside)
         kidlst = node(ffluxptr,mkid)

c        ## debugging output
         if (uprint) then
           write(outunit,101) icrse,jcrse,kcrse,
     .            (val(ivar,icrse,jcrse,kcrse),ivar=1,nvar)
 101       format(" old ",1x,3i4,5e15.7)
         endif

         if (mcapa .gt. 0) then
c            # capacity array:  need to divide by capa in each cell.
c            # modify sgnm which is reset for each grid cell.
c            # Note capa is stored in aux(icrse,jcrse,kcrse,mcapa)
             sgnm = sgnm / alloc(iaddaux(icrse,jcrse,kcrse))
         endif

         do 20 ivar = 1,nvar
            val(ivar,icrse,jcrse,kcrse) = val(ivar,icrse,jcrse,kcrse) +
     1      sgnm*alloc(kidlst+nvar*(lkid-1)+ivar-1)/volume
 20      continue
         iused(icrse,jcrse,kcrse) = iused(icrse,jcrse,kcrse) + norm

c        ## debugging output
         if (uprint) then
           write(outunit,102) mkid,
     .         (val(ivar,icrse,jcrse,kcrse),ivar=1,nvar)
 102       format(" new ","(grid",i3,")",5e15.7)
         endif

 40   continue
c
 99   return
      end
