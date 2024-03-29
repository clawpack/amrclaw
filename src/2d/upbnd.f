c
!> Do conservation fix-up for cells on grid **mptr** that border 
!! finer grids.
!! Do this by adding flux differences stored
!! with each of the fine grids to these cells on grid **mptr**. 
!! 
!! \param[in,out] listbc an array that stores information of all cells that
!! border a finer grids in grid **mptr**. 
!! \param[in] val pointer to the solution on grid **mptr**
!! \param[in] nvar number of equations for the system
!! \param[in] naux number of auxiliary variables
!! \param[in] mitot total number of cells in *i* direction on grid **mptr**
!! \param[in] mjtot total number of cells in *j* direction on grid **mptr**
!! \param[in] maxsp maximum number of segments **lisbc** can describe
!! \param[in] mptr the grid being processed
c ------------------------------------------------------------
c
       subroutine upbnd(listbc,val,nvar,naux,mitot,mjtot,
     1                  maxsp,mptr)
c     1                  maxsp,iused,mptr)
 
      use amr_module
      implicit double precision (a-h,o-z)

 
       dimension val(nvar,mitot,mjtot), iused(mitot,mjtot)
       double precision listbc(5,maxsp)

c  OLD INDEXING
c      iaddaux(i,j) = locaux + i-1 +  mitot*(j-1) 
c    1                 + mitot*mjtot*(mcapa-1)
c  NEW INDEXING - SWITCHED ORDERING
       iaddaux(i,j) = locaux + mcapa-1 +  naux*(i-1) 
     1                 + mitot*naux*(j-1)
 
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
 
      do j=1,mjtot
      do i=1,mitot
         iused(i,j) = 0.
      end do
      end do
 
      locaux = node(storeaux,mptr)
      levc   = node(nestlevel,mptr)
      area   = hxposs(levc)*hyposs(levc)


      if (uprint) write(outunit,*)" upbnding grid ",mptr

      do 40 ispot = 1,maxsp
         icrse = listbc(1,ispot)
         if (icrse.eq.0) go to 99
 
         jcrse = listbc(2,ispot)
         iside = listbc(3,ispot)
c        continue to use iside1/norm for debugging, but should soon delete         
c        this if/then/else block needed due to new categories corresponding
c        to mapped bcs. should still only have one update per side of coarse cell though
         if (iside .lt. 5) then   
           iside1 = iside
         elseif (iside .eq. 5) then
           iside1 = 2
         else  ! iside is 6
           iside1 = 4
         endif
         norm = 2**(iside1-1)
         iflag =iused(icrse,jcrse)/norm
         if (mod(iflag,2).eq.1) then
            write(6,*)" ***  double flux update CAN happen in upbnd ***"
            go to 40
         endif
         mkid = listbc(4,ispot)
         kidlst = node(ffluxptr,mkid)
         lkid = listbc(5,ispot)
c         if (mod(iside,4).gt.1) then
c         modified to include other side options
          if (iside .eq. 2 .or. iside .eq. 3 .or. iside .eq. 6) then
c           (iside .eq. 2 .or. iside .eq. 3)
            sgnm = -1.
         else
c           (iside .eq. 4 .or. iside .eq. 1)
            sgnm = 1.
         endif

c        ## debugging output
         if (uprint) then
           write(outunit,101) icrse,jcrse,
     .         (val(ivar,icrse,jcrse),ivar=1,nvar)
 101       format(" old ",1x,2i4,4e15.7)
         endif

         if (mcapa .gt. 0) then
c            # capacity array:  need to divide by capa in each cell.
c            # modify sgnm which is reset for each grid cell.
c            # Note capa is stored in aux(icrse,jcrse,mcapa)
             sgnm = sgnm / alloc(iaddaux(icrse,jcrse))
         endif

         do 20 ivar = 1,nvar
            val(ivar,icrse,jcrse) = val(ivar,icrse,jcrse) +
     1      sgnm*alloc(kidlst+nvar*(lkid-1)+ivar-1)/area
 20      continue
         iused(icrse,jcrse) = iused(icrse,jcrse) + norm

c        ## debugging output
         if (uprint) then
           write(outunit,102) mkid,
     .         (val(ivar,icrse,jcrse),ivar=1,nvar)
 102       format(" new ","(grid",i3,")",4e15.7)
         endif

 40   continue
c
 99   return
      end
