c
c --------------------------------------------------------------
c
      subroutine stst1
c
      implicit double precision (a-h,o-z)

      include  "call.i"

c
c :::::::::::::::::::::::::::::: STST1 ::::::::::::::::::::::::::::::::
c    intialize a few variables needed before calling user set up
c    routine domain.
c    the spatial and temporal stepsizes are set. the node array
c    is kept as a linked list of free nodes.  "ndfree" points to the
c    head of the list, i.e. first free node.  use first row of each
c    col to hold this pointer, set by the macro "nextfree".
c    the free space list, managed in lfree, will have first and
c    last positions filled with an allocation of zero words,
c    to avoid boundary cases.
c ::::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::::::
c
      ndfree = 1
      do 10 i   = 1, maxgr
         node(nextfree,i) = i+1
 10   continue
c
c the last free node will have a null pointer
 
      node(nextfree, maxgr) = null
c
      lfine = 1
! old style alloc creation
!--      do 20 i  = 1, memsize
!--        alloc(i) = 0.0d0
!-- 20   continue
! new style uses dynamic memory
c
c     Initialize dynamic memory. sets memsize used later in this proc
      call init_alloc()
c
c  initialize linked list of alloc storage as well.
c  first and last locations are dummy placeholders of zero words
c  of allocation each, to avoid boundary cases.
c
      do  40 i  = 1, lfdim
        lfree(i,1) = 0
        lfree(i,2) = 0
 40   continue
      lfree(3,1) =memsize + 2
      lfree(2,1) = 1
      lfree(2,2) =memsize
      lenf       = 3
c
c after kcheck integrations of parent grid, move its refinements.
c finest level grid never needs to have its finer subgrids moved.
c
      do 60 i   = 1, maxlv
         iregridcount(i) = 0
         avenumgrids(i) = 0
         numgrids(i) = 0
         numcells(i) = 0
         lstart(i) = 0
 60      icheck(i) = 0
c
c finish initializing spatial and counting arrays
c
      level      = 2
 70   if (level .gt. mxnest) go to 80
          hxposs(level) = hxposs(level-1) / dfloat(intratx(level-1))
          hyposs(level) = hyposs(level-1) / dfloat(intraty(level-1))
          hzposs(level) = hzposs(level-1) / dfloat(intratz(level-1))
          level         = level + 1
      go to 70
 80   continue

      return
      end
