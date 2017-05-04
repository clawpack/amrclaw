c
c ------------------------------------------------------------
c
       subroutine upbnd(listbc,val,nvar,naux,mitot,
     1                  maxsp,mptr)
c     1                  maxsp,iused,mptr)
 
      use amr_module
      implicit double precision (a-h,o-z)

 
       dimension val(nvar,mitot),listbc(4,maxsp),
     1           iused(mitot)

c  OLD INDEXING
c      iaddaux(i,j) = locaux + i-1 +  mitot*(j-1) 
c    1                 + mitot*mjtot*(mcapa-1)
c  NEW INDEXING - SWITCHED ORDERING
       iaddaux(i) = locaux + mcapa-1 +  naux*(i-1)
 
c
c :::::::::::::::::::::::::::: UPBND :::::::::::::::::::::::::::::
c We now correct the coarse grid with the flux differences stored
c with each of the fine grids. We use an array   iused
c to indicate whether the flux has been updated or not for that zone.
c iused(i) = sum from (l=1,2) i(l)*2**(l-1), where i(l) = 1 if the
c flux for the  l-th side of the (i)-th cell has already been
c updated, and i(l) = 0 if not.
 
c if there is a capacity fn. it needs to be included in update formula
c indicated by mcapa not zero (is index of capacity fn.)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      do 10 i=1,mitot
         iused(i) = 0.
 10   continue
 
      locaux = node(storeaux,mptr)
      levc   = node(nestlevel,mptr)
      area   = hxposs(levc)


      if (uprint) write(outunit,*)" upbnding grid ",mptr

      do 40 ispot = 1,maxsp
         icrse = listbc(1,ispot)
         if (icrse.eq.0) go to 99

         iside = listbc(2,ispot)
c        continue to use iside1/norm for debugging, but should soon delete         
c        this if/then/else block needed due to new categories corresponding
c        to mapped bcs. should still only have one update per side of coarse cell though
         if (iside .lt. 3) then
           iside1 = iside
         else  ! should not happen
           write(*,*) "*** Error in upbud.f ***"
         endif
         norm = 2**(iside1-1)
         iflag =iused(icrse)/norm
         if (mod(iflag,2).eq.1) then
            write(6,*)" ***  double flux update CAN happen in upbnd ***"
            go to 40
         endif
         mkid = listbc(3,ispot)
         kidlst = node(ffluxptr,mkid)
         lkid = listbc(4,ispot)
c         modified to include other side options
          if (iside .eq. 2) then
c           (iside .eq. 2)
            sgnm = -1.
         else
c           (iside .eq. 1)
            sgnm = 1.
         endif

c        ## debugging output
         if (uprint) then
           write(outunit,101) icrse,
     .         (val(ivar,icrse),ivar=1,nvar)
 101       format(" old ",1x,i4,4e15.7)
         endif

         if (mcapa .gt. 0) then
c            # capacity array:  need to divide by capa in each cell.
c            # modify sgnm which is reset for each grid cell.
c            # Note capa is stored in aux(icrse,mcapa)
             sgnm = sgnm / alloc(iaddaux(icrse))
         endif

         do 20 ivar = 1,nvar
            val(ivar,icrse) = val(ivar,icrse) +
     1       sgnm*alloc(kidlst+nvar*(lkid-1)+ivar-1)/area
 20      continue
         iused(icrse) = iused(icrse) + norm

c        ## debugging output
         if (uprint) then
           write(outunit,102) mkid,
     .         (val(ivar,icrse),ivar=1,nvar)
 102       format(" new ","(grid",i3,")",4e15.7)
         endif

 40   continue
c
 99   return
      end
