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
      bndListSize = 8*maxgr   ! guess at initial size

      ! allocate the bndList
      if (.not. allocated(bndList)) then
         allocate(bndList(bndListSize,2))
         print *, "bndList allocated..."
      else
         print *, "bndList already allocated!"
      endif

      ! thread the bndry list
      do i = 1, bndListSize
         bndList(i,nextfree) = i+1
      end do

      bndList(bndListSize,nextfree) = null
      ndfree_bnd = 1

      end
c
c -----------------------------------------------------------------
c
!> Preprocess each grid on level **level** to have a linked list of 
!! other grids at the same level that supply ghost cells.
!!
!! The linked list is implemented with an array.
!! Each node in this linked list is represented by a row (with 2 elements) 
!! in the array, named **bndList**.
!!
!! bndList(pos,gridNbor) is grid number stored at node **pos**
!! bndList(pos,nextfree) is pointer to next node in the linked list
!!
!! node(bndListSt,mptr) point to the head node in this linked list.
!! Thus bndList(node(bndListSt,mptr),gridNbor) is grid number 
!! stored at first node of the linked list.
!!
!!
      subroutine makeBndryList(level)
c
      use amr_module
      implicit none

      integer level, n, levSt, k, nborCount
      integer nodget_bnd, nextSpot, prevNbor, msrc, mptr
      integer imin, imax, jmin, jmax
      integer imlo, imhi, jmlo, jmhi
      integer ixlo, ixhi, jxlo, jxhi

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
         jmin = node(ndjlo,mptr) - nghost  ! may also use for filval stuff,
         jmax = node(ndjhi,mptr) + nghost  ! change nghost to mbuff, etc
         nborCount = 0
         
         do k = 1, numgrids(level)  ! loop over all other grids once to find touching ones 
            if (k .eq. n) cycle     ! dont count yourself as source grid
            msrc = listOfgrids(levSt+k-1)

            ! Check if grid mptr and patch intersect
            imlo = node(ndilo, msrc)
            jmlo = node(ndjlo, msrc)
            imhi = node(ndihi, msrc)
            jmhi = node(ndjhi, msrc)

            ixlo = max(imlo,imin)
            ixhi = min(imhi,imax)
            jxlo = max(jmlo,jmin)
            jxhi = min(jmhi,jmax)

            if (ixlo .le. ixhi .and. jxlo .le. jxhi) then ! put on bnd list for mptr
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
!> Free the linked list of intersecting "boundary" grids for grid 'mold'
!! that is no longer active.
!! The linked list starts at node(bndListSt, mold).
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
