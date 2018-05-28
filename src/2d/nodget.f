c
c ------------------------------------------------------------
c
!> Get first free node of the linked list kept in node
!! array. adjust pointers accordingly.
!!
!! This function is used to create a new grid descriptor, like
!! mptr = new grid_class in c++.
!!
      integer function nodget()
c
      use amr_module
      implicit double precision (a-h,o-z)
      integer maxgrIncrement/10000/

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
 101         format("    level ",i4," has ",i7,' grids')
          end do
          write(*,*)" Could need twice as many grids as on any given"
          write(*,*)" level if regridding/birecting"
!          stop
!         new way gets more storage and continues
          istatus = 1
          call resize_nodes(maxgr + maxgrIncrement, istatus)
          if (istatus > 0) then
             write(*,102) maxgr,maxgrIncrement
102          format(' unable to increase nodal space by ',i8,' grids',/,
     .              'resize failed')
             stop
          endif
c
c  allocate next node, update pointers
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

          !stop
          call resize_bndryList()
c
c     ## adjust pointers for next bndry entry
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
!> Make (modify) array of grid numbers, listOfGrids, (after sorting them 
!! in the linked list so they are in decreasing order of workload, 
!! done in arrangeGrid()).
!!
!! This is done every time there is regridding, initial gridding,
!! or restarting.  Most often finest level is regridded, so
!! put it last in array. **lbase** is the level that didnt change, so 
!! only redo from lbase+1 to lfine.
!! \param[in] lbase all levels from **lbase**+1 to the finest get
!! modifed in the array, listOfGrids
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
