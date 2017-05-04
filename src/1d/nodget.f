c
c ------------------------------------------------------------
c
      integer function nodget()
c
      use amr_module
      implicit double precision (a-h,o-z)

c
c ::::::::::::::::: NODGET ::::::::::::::::::::::::::::::::::::;
c nodget =  get first free node of the linked list kept in node
c            array. adjust pointers accordingly.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      if (ndfree .ne. null) go to 10
          write(outunit,100) maxgr
          write(*,100)       maxgr
100       format(' out of nodal space - allowed ',i8,' grids')
          do level = 1, lfine
             write(*,101) level,numgrids(level)
 101         format("    level ",i4," has ",i6,'grids')
          end do
          write(*,*)" Could need twice as many grids as on any given"
          write(*,*)" level if regridding/birecting"
          stop
c
c  update pointers
c
 10     nodget         = ndfree
        ndfree         = node(nextfree,ndfree)
c
c  initialize new  block
c
        do 20 i        = 1, nsize
           node(i,nodget) = 0
 20     continue
c
        do 30 i         = 1, rsize
           rnode(i,nodget) = 0.0d0
 30     continue
c
      return
      end
c
c ------------------------------------------------------------
c
      integer function nodget_bnd()
c
      use amr_module
      implicit double precision (a-h,o-z)

c
c ::::::::::::::::: NODGET_BND ::::::::::::::::::::::::::::::::::::;
c nodget_bnd =  same as above but for bndry list
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
c
      if (ndfree_bnd .ne. null) go to 10
          write(outunit,100) bndListSize
          write(*,100)       bndListSize
100       format(' out of bndry space - allowed ',i5,' bndry grids')
          ! calc average number of bndry nbors per grid
          nborTotal = 0
          numGridsTotal = 0
          do lev = 1, lfine
            numGridsTotal = numGridsTotal + numgrids(lev)
            do mptr = 1, numgrids(lev)
               nborTotal = nborTotal + node(bndListNum,mptr)  
            end do
          end do
          avgNbors = float(nborTotal)/numgridsTotal
          write(*,101) numGridsTotal,nborTotal,avgNbors
 101      format(" There are ",i8," total grids", i10," bndry nbors",
     .           " average num/grid ",f10.3)

          stop
c
c     ## adjust pointers
c
 10   nodget_bnd      = ndfree_bnd
      ndfree_bnd      = bndList(ndfree_bnd,nextfree)
c
c     ##  initialize to 0
c
      bndList(nodget_bnd,1) = 0
      bndList(nodget_bnd,2) = 0
c
      return
      end
c
c -----------------------------------------------------------------
c
      subroutine makeGridList(lbase)
c
      use amr_module
      implicit none

      integer lbase, levSt, lev, mptr, n

c :::::::::::::::::::::::::::: make_gridList :::::::::::::::::::::::::
c     make array of grid numbers (after sorting them so in decreasing
c     order of workload, done in arrangeGrid and put back into linked 
c     list. Done every time there is regridding, initial gridding,
c     or restarting.  Most often finest level is regridded, so
c     put it last in array. lbase is the level that didnt change, so 
c     only redo from lbase+1 to lfine.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      !write(*,*)"mgl: lbase,lfine",lbase,lfine
      do lev = lbase+1, lfine
         levSt = listStart(lev) 
         !write(*,*)"mgl: level ",lev," starts at ",levSt
         mptr = lstart(lev)
c        traverse linked list into array. list already sorted by arrangegrids
         do n = 1, numgrids(lev)
            listOfGrids(levSt+n-1) = mptr
            mptr = node(levelptr,mptr)
         end do
c
c        next level starts one after where this one ends.
c        Using a sentinel in dimension of
c        listStart so no need to test if level = mxnest
c
         listStart(lev+1) = levSt + numgrids(lev)
      end do

      return
      end
c
c -----------------------------------------------------------------
c
      subroutine initBndryList()

      use amr_module
      implicit none
      
      integer i
c
c  need to manage the boundary List too
c
      do i = 1, bndListSize
         bndList(i,nextfree) = i+1
      end do

      bndList(bndListSize,nextfree) = null
      ndfree_bnd = 1

      end
c
c -----------------------------------------------------------------
c
      subroutine makeBndryList(level)
c
      use amr_module
      implicit none

      integer level, n, levSt, k, nborCount
      integer nodget_bnd, nextSpot, prevNbor, msrc, mptr
      integer imin, imax
      integer imlo, imhi
      integer ixlo, ixhi

c :::::::::::::::::::::::::::: makeBndryList :::::::::::::::::::::::::
c     preprocess each grid to have linked list of other grids at
c     same level that supply ghost cells.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     traverse linked list into array. list already sorted by arrangegrids
      levSt = listStart(level) 
      do n = 1, numgrids(level)
         mptr = listOfGrids(levSt+n-1)
         imin = node(ndilo,mptr) - nghost  ! ghost cells included since
         imax = node(ndihi,mptr) + nghost  ! this is what you want to fill
                                           ! may also use for filval stuff,
                                           ! change nghost to mbuff, etc
         nborCount = 0
         
         do k = 1, numgrids(level)  ! loop over all other grids once to find touching ones 
            if (k .eq. n) cycle     ! dont count yourself as source grid
            msrc = listOfgrids(levSt+k-1)

            ! Check if grid mptr and patch intersect
            imlo = node(ndilo, msrc)
            imhi = node(ndihi, msrc)

            ixlo = max(imlo,imin)
            ixhi = min(imhi,imax)

            if (ixlo .le. ixhi) then ! put on bnd list for mptr
               nborCount = nborCount + 1
               nextSpot = nodget_bnd()   
               bndList(nextSpot,gridNbor) = msrc
               ! get spot in bnd list. insert next grid at front to avoid traversing
               bndList(nextSpot,nextfree) =  node(bndListSt,mptr)
               node(bndListSt,mptr) = nextSpot
            endif

         end do

!        save final count
         node(bndListNum,mptr) = nborcount
      end do

      return
      end
c
c -----------------------------------------------------------------
c
      subroutine freeBndryList(mold)
c
      use amr_module
      implicit none

      integer nborCount, mold,nextSpot, i, nextnext

c :::::::::::::::::::::::::::: freeBndryList :::::::::::::::::::::::::
c     free the linked list of intersecting "boundary" grids for grid 'mold'
c     that is no longer active
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

           nborCount = node(bndListNum,mold) ! count for this grid
           nextSpot  = node(bndListSt,mold)  ! first index of this grids nbors
           do i = 1, nborCount
               nextnext = bndList(nextSpot,nextfree)
               call putnod_bnd(nextSpot)
               nextSpot = nextnext
           end do

      return
      end
