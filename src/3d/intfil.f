c
c ------------------------------------------------------
c
       subroutine intfil(val,mi,mj,mk,time,flaguse,nrowst,ncolst,nfilst,
     2                   ilo,ihi,jlo,jhi,klo,khi,level,nvar,naux,msrc)
c
c ::::::::::::::::::::::: INTFIL ::::::::::::::::::::::::::::::::;
c  INTFIL: interpolates values for a patch at the specified level and
c  location, using values from grids at LEVEL and coarser, if nec.
c
c  take the intersection of a grid patch with corners at
c  ilo,ihi,jlo,jhi,klo,khi
c  and all grids mptr at LEVEL.  If there is a non-null intersection
c  copy the solution vaues from mptr (at TIME) into VAL array.
c  assumes patch at same level so do straight copy, not skipping
c  every intrat or doing any interpolation here,
c  assume called in correct order of levels, so that when copying
c  is ok to overwrite.
c
c  N.B.: there are no dummy points around patch, since
c        this is not an official "grid" but a "patch".
c
c  used array marks when point filled. filpatch checks if points left over 
c  after intersections at specified level.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      use amr_module
      implicit double precision (a-h,o-z)

c
      dimension   val(nvar,mi,mj,mk)
      logical     tinterp
c
      iadd(ivar,i,j,k)   = loc    +    (ivar-1)
     &                            +    (i-1)*nvar
     &                            +    (j-1)*nvar*mitot
     &                            +    (k-1)*nvar*mitot*mjtot
      iadnew(ivar,i,j,k) = locnew +    (ivar-1)
     &                            +    (i-1)*nvar
     &                            +    (j-1)*nvar*mitot
     &                            +    (k-1)*nvar*mitot*mjtot
      iadold(ivar,i,j,k) = locold +    (ivar-1)
     &                            +    (i-1)*nvar
     &                            +    (j-1)*nvar*mitot
     &                            +    (k-1)*nvar*mitot*mjtot
      dimension flaguse(ilo:ihi, jlo:jhi, klo:khi)
c
      dt     = possk(level)
      teps   = dt / 10.

      nrow = ihi - ilo + 1
      ncol = jhi - jlo + 1
      nfil = khi - klo + 1
      do k = klo, khi
         do j = jlo, jhi
            do i = ilo, ihi
               flaguse(i,j,k) = 0.0
            end do
         end do
      end do

!      mptr   = lstart(level)
! 10   if (mptr .eq. 0) go to 105
       
      if (msrc .eq. -1) then ! use array of all grids at this level
         numg = numgrids(level)
         levSt = listStart(level)
      else  ! use new bndry list for each grid if at same level
         bndLoc = node(bndListSt,msrc)  ! index of first grid in bndList
         bndNum = node(bndListNum,msrc)
         nextSpot = node(bndListSt, msrc) ! initialize
         numg = bndNum
      endif

      do icount = 1, numg   ! combining old and new way to loop over grids
      
         if (msrc .eq. -1) then
            mptr = listOfGrids(levSt+icount-1)
         else
            mptr = bndList(nextSpot,gridNbor)
         endif
c
c     :::  check if grid mptr and patch intersect
c
      imlo = node(ndilo, mptr)
      jmlo = node(ndjlo, mptr)
      kmlo = node(ndklo, mptr)
      imhi = node(ndihi, mptr)
      jmhi = node(ndjhi, mptr)
      kmhi = node(ndkhi, mptr)

      nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
      nz     = node(ndkhi,mptr) - node(ndklo,mptr) + 1
      mitot = nx + 2*nghost
      mjtot = ny + 2*nghost
      mktot = nz + 2*nghost

      ixlo = max(imlo,ilo)
      ixhi = min(imhi,ihi)
      jxlo = max(jmlo,jlo)
      jxhi = min(jmhi,jhi)
      kxlo = max(kmlo,klo)
      kxhi = min(kmhi,khi)
c
      if (.not.((ixlo .le. ixhi .and. jxlo .le. jxhi) .and.
     &          (kxlo .le. kxhi))) go to 90
c
c  :::  grids intersect. figure out what time to use.
c  :::  alphai = 1 for new time; 0 for old time
c
c     alphai = (time+dt - rnode(timemult,mptr)) / dt
c     alphac = 1.d0-alphai

      alphac = (rnode(timemult,mptr)-time)/dt
      alphai = 1.d0-alphac

      if ((alphai .lt. -teps) .or. (alphai .gt. 1.d0+teps)) then
          write(outunit,900) time, mptr, level
          write(*,900) time, mptr, level
 900      format(' time wanted ',e15.7,' not available from grid ',i4,
     1           'level',i4)
          write(outunit,901) ilo,ihi,jlo,jhi,klo,khi,mptr,level,time,
     .                 rnode(timemult,mptr),alphai,teps
          write(*,901) ilo,ihi,jlo,jhi,klo,khi,mptr,level,time,
     .                 rnode(timemult,mptr),alphai,teps
          call outtre(mstart,.false.,nvar,naux)
          stop
      endif
c
      tinterp = .false.
      if (dabs(alphai - 1.d0) .lt. teps) then
          loc = node(store1,mptr)
      else if (dabs(alphai) .lt. teps) then
          loc = node(store2,mptr)
          if (level .eq. mxnest) then
             write(outunit,901) ilo,ihi,jlo,jhi,klo,khi,mptr,level,time,
     .                          rnode(timemult,mptr),alphai,teps
             write(*,901) ilo,ihi,jlo,jhi,klo,khi,mptr,level,time,
     .                          rnode(timemult,mptr),alphai,teps
             stop
           endif
      else
          locold  = node(store2,mptr)
          locnew  = node(store1,mptr)
          tinterp = .true.
          if (level .eq. mxnest) then
             write(outunit,901) ilo,ihi,jlo,jhi,klo,khi,mptr,level,time,
     .                    rnode(timemult,mptr),alphai,teps
             write(*,901) ilo,ihi,jlo,jhi,klo,khi,mptr,level,time,
     .                    rnode(timemult,mptr),alphai,teps
             stop
          endif
      endif
 901  format(' trying to interpolate from previous time values ',/,
     .       ' for a patch with corners ilo,ihi,jlo,jhi,klo,khi:'
     .       ,/,2x,6i10,/,
     .       ' from source grid ',i4,' at level ',i4,/,
     .       ' time wanted ',e15.7,' source time is ',e15.7,/,
     .       ' alphai, teps ',2e15.7)
c
      if (.not. tinterp) then
c     ::: no time interp. copy the solution values
         do 35 k = kxlo, kxhi
         do 25 j = jxlo, jxhi
         do 15 i = ixlo, ixhi
         do 45 ivar = 1, nvar
             val(ivar,i-ilo+nrowst,j-jlo+ncolst,k-klo+nfilst) =
     1            alloc(iadd(ivar,i-imlo+nghost+1,
     2                       j-jmlo+nghost+1,
     3                       k-kmlo+nghost+1))
 45          continue
             flaguse(i,j,k) = 1.d0
 15      continue
 25      continue
 35      continue
      else
c     ::: linear interpolation in time 
         do 75 k = kxlo, kxhi
         do 65 j = jxlo, jxhi
         do 55 i = ixlo, ixhi
         do 85 ivar = 1, nvar
             val(ivar,i-ilo+nrowst,j-jlo+ncolst,k-klo+nfilst) =
     1            alloc(iadnew(ivar,i-imlo+nghost+1,
     2                         j-jmlo+nghost+1,
     3                         k-kmlo+nghost+1))*alphai +
     1            alloc(iadold(ivar,i-imlo+nghost+1,
     2                         j-jmlo+nghost+1,
     3                         k-kmlo+nghost+1))*alphac
 85         continue
            flaguse(i,j,k) = 1.d0
 55      continue
 65      continue
 75      continue
      endif
c
!90   mptr = node(levelptr, mptr)
 90      continue

         if (msrc .ne. -1) then
            nextSpot = bndList(nextSpot,nextfree)
         endif
      end do
!      go to 10
c
 105  continue
 
c  set used array points which intersect boundary to be equal to 1;
c  they will be set elsewhere
 
      if (khi .ge. kregsz(level)) then
        do 1100 k = max(kregsz(level),klo), khi
        do 1100 j = jlo, jhi
        do 1100 i = ilo, ihi
           flaguse(i,j,k) = 1.d0
1100    continue
      endif
 
      if (klo .lt. 0) then
        nfilend = nfilst + nfil - 1
        do 1200 k = klo, min(-1,nfilend)
        do 1200 j = jlo, jhi
        do 1200 i = ilo, ihi
            flaguse(i,j,k) = 1.d0
1200    continue
      endif

      if (jhi .ge. jregsz(level)) then
        do 1300 k = klo, khi
        do 1300 j = max(jregsz(level),jlo), jhi
        do 1300 i = ilo, ihi
           flaguse(i,j,k) = 1.d0
1300    continue
      endif
 
      if (jlo .lt. 0) then
        ncolend = ncolst + ncol - 1
        do 1400 k = klo, khi
        do 1400 j = jlo, min(-1,ncolend)
        do 1400 i = ilo, ihi
           flaguse(i,j,k) = 1.d0
1400    continue
      endif
 
      if (ilo .lt. 0) then
        nrowend = nrowst + nrow - 1
        do 1500 k = klo, khi
        do 1500 j = jlo, jhi
        do 1500 i = ilo, min(-1,nrowend)
           flaguse(i,j,k) = 1.d0
1500    continue
      endif
 
      if (ihi .ge. iregsz(level)) then
        do 1600 k = klo, khi
        do 1600 j = jlo, jhi
        do 1600 i = max(iregsz(level),ilo), ihi
           flaguse(i,j,k) = 1.d0
1600    continue
      endif
c
      return
      end
