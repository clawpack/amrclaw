c
c -----------------------------------------------------------
c
      subroutine colate2 (badpts, len, lcheck, nUniquePts, lbase)
c
      use amr_module
      implicit  double precision (a-h,o-z)
      dimension badpts(2,len)
      dimension ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)

c
c    index for flag array now based on integer index space, not 1:mibuff,1:mjbuff
c    but if grid extends outside domain, not supposed to have flagged points
      iadd(i,j) = locamrflags + i-(ilo-mbuff) + mibuff*(j-(jlo-mbuff))
c
c
c *************************************************************
c
c colate2 = takes each grids flagged points at level lcheck
c          and puts their (i,j) cell centered
c          indices into the badpts array.
c          To insure proper nesting, must get rid of flagged points
c          that dont fit into properly nested domain. Grids
c          with flagged points include buffered region (size mbuff)now.
c          THIS NEW VERSION may have duplicate points. need to sort
c          and remove when colating.
c
c if checking flagged pt for nesting is expensive, might consider not doing it
c and revising projec2 instead. if new fine grid not nested, then go through
c flagged points and throw out. But presumably many grids will make it through
c without having to check all points.
c
c *************************************************************
c
c  any cell strictly interior (one inside boundary) to the grid proper  is already properly nested.
c so if it is flagged, even if parent grid is at level lbase (not moving), it is ok.
c Potential problem is only from buffered cells exterior of grid, in the mbuff
c region, and the first interior row/col of cells. So only need to check those.
c if base level is level 1, then must be ok, for convex domain.
c
       mbuff = max(nghost,ibuff)  ! new way of expanding grid to do buffering in place
       index = 0  ! for putting into badpts array
       lratiox = 1
       lratioy = 1
       do lev = lbase,lcheck-1  
          lratiox = lratiox * intratx(lev)
          lratioy = lratioy * intraty(lev)
       end do


      mptr = lstart(lcheck)
 10      if (node(numflags,mptr) .eq. 0) go to 70    !simple bypass if no tags
         write(outunit,*)" colating flags on grid ",mptr

c        handle each of 4 sides (in 2D)
c        set tags to negative val. reset to positive if they have a home     
         ilo = node(ndilo,mptr)
         ihi = node(ndihi,mptr)
         jlo = node(ndjlo,mptr)
         jhi = node(ndjhi,mptr)
         nx = ihi - ilo + 1
         ny = jhi - jlo + 1
         mibuff = nx + 2 *mbuff
         mjbuff = ny + 2 *mbuff


         locamrflags = node(storeflags,mptr)
         call setNeg(alloc(locamrflags),ilo,ihi,jlo,jhi,mbuff,ico)
         if (ico .eq. 0) go to 70  ! no points flagged on this grid; nothing to do

c        take care of removing flagged pts at solid bndry in buffer zone. those pts not nested
c        (but it doesnt matter). The last row of flagged pts on grid touching bndry is considered
c        nested, since touches domain bndry.
         call setPhysBndry(alloc(locamrflags),ilo,ihi,jlo,jhi,
     .                     mbuff,lcheck)

c        next: deal with periodicity by looping over 3 regions in each coordinate dir
c        for each non-empty region for grid mptr (enlarged by mbuff) need to intersect with base grids  
c        if periodic domain buffered region can extend out
c        if not, reset to  physical region only

         call  setIndices(ist,iend,jst,jend,ilo-mbuff,ihi+mbuff,
     &                   jlo-mbuff,jhi+mbuff,ishift,jshift,lcheck)

         if (xperdom) then  ! no need to check against phys bndry
             ilo_use = ilo - mbuff
             ihi_use = ihi + mbuff  
         else  !  isnt there a way to put this in setIndices?
            ilo_use = max(ilo-mbuff,0)
            ihi_use = min(ihi+mbuff,iregsz(lcheck)-1)
         endif
         if (yperdom) then  ! no need to check against phys bndry
             jlo_use = jlo - mbuff
             jhi_use = jhi + mbuff  
         else
            jlo_use = max(jlo-mbuff,0)
            jhi_use = min(jhi+mbuff,jregsz(lcheck)-1)
         endif

         do 66 ireg = 1, 3
            i1 = max(ilo_use, ist(ireg))
            i2 = min(ihi_use, iend(ireg))
         do 65 jreg = 1, 3
            j1 = max(jlo_use, jst(jreg))
            j2 = min(jhi_use, jend(jreg))           
            if (.not.((i1 .le. i2) .and. (j1 .le. j2))) go to 65  ! no patch in this pregion
c
c           otherwise need to intersect wrapped patch with base grids 

            mbase = lstart(lbase)
 40         continue
          
c           search neighboring grids to see if bndry slabs nested
c           scale up their coords then compare to current grid         
c           add 1 to right side to anchor properly. 
c           then add/sub 1 to shrink, for proper nesting
c           NB: this alg. will throw out a point at edge of nested grid since it
c           doesnt check for adjacent grid. Better soln (for the future) to
c           create domflags for each grid, (not whole domain)
c
               iblo =  node(ndilo, mbase)*lratiox + 1
               ibhi = (node(ndihi, mbase)+1)*lratiox - 2
               jblo =  node(ndjlo, mbase)*lratioy + 1
               jbhi = (node(ndjhi, mbase)+1)*lratioy - 2
c
c NB: intersecting base grid (shrunken by 1 in this new alg. as way to make sure of nesting)
c     with enlarged mptr grid - to see if flagged cells in buffer region are really ok after all
            ixlo = max(iblo,i1+ishift(ireg))   ! NB using periodically shifted coords, if nec.
            ixhi = min(ibhi,i2+ishift(ireg))
            jxlo = max(jblo,j1+jshift(jreg))
            jxhi = min(jbhi,j2+jshift(jreg))            
            if (.not.((ixlo .le. ixhi) .and. (jxlo .le. jxhi))) go to 50
c
c           grid intersects, set flags to abs val. 
c           ### easier to loop over whole grid, but should write sub to to bdry vals only
            do 45 j = jxlo-jshift(jreg), jxhi-jshift(jreg)
            do 45 i = ixlo-ishift(ireg), ixhi-ishift(ireg)
!--              if (i .ge. ilo .and. i .le. ihi .and. j .ge. jlo
!--     &             .and. j .le.  jhi) go to 43  ! no printing for debugging purposes
!--              write(outunit,949) i,j,mptr
!-- 949          format(" setting pt. ",2i5," in grid ",i5," to positive")
 43           alloc(iadd(i,j)) = abs(alloc(iadd(i,j)))
 45         continue

 50         mbase = node(levelptr, mbase)
            if (mbase .ne. 0) go to 40

c
c           insert (wrapped) flagged points into badpts array
c           so dont have to recompute shifted indices in another loop
c           put remaining (positive) flags into badpts array. 
c           at this point still has duplicates
c           give points the indices from integer region space.
           do 60 j   = j1, j2
           do 60 i   = i1, i2
             if (alloc(iadd(i,j)) .gt. goodpt) then  ! neg means no home was found. throw out
              index = index + 1
c  WARNING: to match orig program note we ADD .5, not subtract. old program used 1 based indexing
c  for grid flagging array. we are using 0 based, so need to add to match previous
c  grid fitting (dont want to change all routines downstream)
              badpts(1,index) = dble(i+ishift(ireg))+.5   ! in case periodic, put flagged buffer pt
              badpts(2,index) = dble(j+jshift(jreg))+.5   ! in badpts in wrapped coords
c             but note, need to loop over coords of enlarged buffered grid, not wrapped grid
c             just inserting points in wrapped coords
c             write(outunit,101) badpts(1,index),badpts(2,index)
             else if (alloc(iadd(i,j)) .lt. goodpt) then
              write(outunit,939) i,j
 939         format("NOT NESTED: ignoring point ",2i5)
 101         format(2f6.1)
             endif
 60        continue

 65         continue
 66         continue



c
 70     continue

c  done colating - safe to reclam
        call reclam(locamrflags,mibuff*mjbuff)
c
        mptr = node(levelptr, mptr)
       if (mptr .ne. 0) go to 10


      npts = index 
      if (gprint) then
        write(outunit,100) npts, lcheck
 100    format( i5,' flagged points initially colated on level ',i4)
      endif
c
c colate flagged points into single integer array for quicksorting
c
      call driveSort(npts,badpts,lcheck,nUniquePts,mbuff)
   

 99   return
      end

c
c -------------------------------------------------------------
c
       subroutine driveSort(npts,badpts,level,index,mbuff)

      use amr_module
      implicit  double precision (a-h,o-z)
      dimension badpts(2,npts)
      dimension iflags(npts), ixArray(npts)
 
      iadd(i,j) = (i+mbuff)  + (isize+2*mbuff)*(j+mbuff)
c
c  convert using one dimensional ordering of badpts array as if
c  it covered entire domain (xprob by yprob) on this level
c
      isize = iregsz(level)
      jsize = jregsz(level)

      do k = 1, npts
        i = badpts(1,k)-.5  ! remember was shifted when put  into array
        j = badpts(2,k)-.5  
        intEquiv = iadd(i,j)
c        write(*,*)i,j," has equivalent integer ",intEquiv
        iflags(k) = intEquiv
      end do

      call qsorti(ixArray, npts, iflags)

c copy back to badpts, in sorted order, removing duplicates      
      k = 1
      index = 0
      do while (k .le. npts) 
         intEquiv = iflags(ixArray(k))
         index = index + 1
         badpts(2,index) = intEquiv/(isize+2*mbuff) + .5 -mbuff
         badpts(1,index) = mod(intEquiv,(isize+2*mbuff)) + .5 -mbuff
c        write(*,101) index, badpts(1,index),badpts(2,index)
c 101    format(" index ",i5," restoring badpts ",2f6.1)
         k = k + 1
         do while ( k.le. npts)    ! skip over duplicates
            if (iflags(ixArray(k)) .eq. iflags(ixArray(k-1))) then
c           write(*,*)" duplicate in sorted array loc ",ixarray(k)
               k = k+1
            else
               exit   ! back to outer loop
            endif
         end do
       if (k .gt. npts) exit !did we stop because we ran off end or pts not equal
      end do

      if (gprint) then
          write(outunit,929) index
 929      format(i5," flagged pts after removing duplicates and ",
     &           " non-nested flags")
      endif

      return
      end
