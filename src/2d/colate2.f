c
c -----------------------------------------------------------
c
      subroutine colate2 (badpts, len, lcheck, nUniquePts)
c
      use amr_module
      implicit  double precision (a-h,o-z)
      dimension badpts(2,len)

      iadd(i,j) = locamrflags + i-ilo-mbuff + mibuff*(j-jlo-mbuff)
c
c
c *************************************************************
c
c colate2 = takes eachs grid flagged points at level lcheck
c          and puts their (i,j) cell centered
c          indices into the badpts array.
c          To insure proper nesting,  get rid of flagged point
c          that dont fit into properly nested domain. Grids
c          with flagged points include buffered region now.
c          THIS NEW VERSION may have duplicate points. need to sort
c          and remove when colating.
c
c *************************************************************
c
c by definition any cell interior to the grid proper is already properly nested.
c so if it is flagged, even if parent grid is at level lbase (not moving), it is ok.
c Potential problem is only from buffered cells exterior of grid, in the mbuff
c region. So only need to check those.
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

       
c        STILL NEED TO HANDLE PERIODIC CASE
         call setPhysBndry(alloc(locamrflags),ilo,ihi,jlo,jhi,
     .                     mbuff,lcheck)

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

            ixlo = max(iblo,ilo)
            ixhi = min(ibhi,ihi)
            jxlo = max(jblo,jlo)
            jxhi = min(jbhi,jhi)
            
            if (.not.((ixlo .le. ixhi) .and. (jxlo .le. jxhi))) go to 50
c
c           grid intersects, set flags to abs val. easier to loop over whole grid
            do 45 i = ixlo, ixhi
            do 45 j = jxlo, jxhi
              alloc(iadd(i,j)) = abs(alloc(iadd(i,j)))
 45         continue

 50         mbase = node(levelptr, mbase)
            if (mbase .ne. 0) go to 40

c       put remaining (positive) flags into badpts array. 
c       at this point still has duplicates
c       give points the indices from integer region space.
        do 20 j   = jlo-mbuff, jhi+mbuff
        do 20 i   = ilo-mbuff, ihi+mbuff
c          if (alloc(iadd(i,j)) .ne. goodpt) then
          if (alloc(iadd(i,j)) .gt. 0.) then  ! neg means no home was found. throw out
            index = index + 1
            badpts(1,index) = dble(i)-.5
            badpts(2,index) = dble(j)-.5
          endif
 20     continue
c
 70   mptr = node(levelptr, mptr)
      if (mptr .ne. 0) go to 10


      npts = index 
      if (gprint) then
        write(outunit,100) npts, lcheck
 100    format( i5,' flagged points initially colated on level ',i4)
      endif
c
c colate flagged points into single integer array for quicksorting
c
      call driveSort(npts,badpts,lcheck,nUniquePts)
   

 99   return
      end

c
c -------------------------------------------------------------
c
       subroutine driveSort(npts,badpts,iflags,level,index)

      use amr_module
      implicit  double precision (a-h,o-z)
      dimension badpts(2,npts)
      dimension iflags(npts), ixArray(npts)
 
      iadd(i,j) = i + isize*j
c
c  convert using one dimensional ordering of badpts array as if
c  it covered entire domain (xprob by yprob) on this level
c
      isize = iregsz(level)
      jsize = jregsz(level)

      do k = 1, npts
        i = badpts(1,k)+.5  ! remember was shifted when put  into array
        j = badpts(2,k)+.5  
        intEquiv = iadd(i,j)
        write(*,*)i,j," has equivalent integer ",intEquiv
        iflags(k) = intEquiv
      end do

      call qsorti(ixArray, npts, iflags)

c copy back to badpts, in sorted order, removing duplicates      
      k = 1
      index = 0
      do while (k .le. npts) 
         intEquiv = iflags(ixArray(k))
         index = index + 1
         badpts(2,index) = intEquiv/isize - .5
         badpts(1,index) = intEquiv mod isize - .5
         do while ( k<npts)    ! skip over duplicates
            if (iflags(ixArray(k)) .eq. iflags(ixArray(k+1))) k = k+1
         end do
       if (k .eq. npts) exit
      end do

      if (gprint) write(*,*) index," flagged pts after",
     .                      " removing duplicates"
      return
      end
