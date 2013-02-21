c
c -------------------------------------------------------------
c
      subroutine spest2 (nvar,naux,lcheck,start_time)
c
      use amr_module
      implicit double precision (a-h,o-z)

      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/
      integer listgrids(numgrids(lcheck))
 
c :::::::::::::::::::::::::: SPEST2 :::::::::::::::::::::::::::::::::::
c For all grids at level lcheck:
c   Call user-supplied routine flag2refine to flag any points where
c   refinement is desired based on user's criterion.  
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
!$    maxthreads = omp_get_max_threads()
      call prepgrids(listgrids,numgrids(lcheck),lcheck)

c       mptr = lstart(lcheck)
c 5     continue
!$OMP PARALLEL DO PRIVATE(jg,mptr,nx,ny,mitot,mjtot,locnew,locaux),
!$OMP&            PRIVATE(time,dx,dy,xleft,ybot,xlow,ylow,locbig),
!$OMP&            PRIVATE(locold,mbuff,mibuff,mjbuff,locamrflags,i),
!$OMP&            SHARED(numgrids,listgrids,lcheck,nghost,nvar,naux),
!$OMP&            SHARED(start_time,possk,flag_gradient,ibuff),
!$OMP&            SHARED(tolsp,alloc,node,rnode,hxposs,hyposs),
!$OMP&            DEFAULT(none),
!$OMP&            SCHEDULE(DYNAMIC,1)
       do  jg = 1, numgrids(lcheck)
          mptr = listgrids(jg)
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          locnew = node(store1,mptr)
          locaux = node(storeaux,mptr)
          time   = rnode(timemult,mptr)
          dx     = hxposs(lcheck)
          dy     = hyposs(lcheck)
          xleft  = rnode(cornxlo,mptr)
          ybot   = rnode(cornylo,mptr)
          xlow   = xleft - nghost*dx
          ylow   = ybot - nghost*dy
c
          locbig = igetsp(mitot*mjtot*nvar)
          node(tempptr,mptr) = locbig
c         # straight copy into scratch array so don't mess up latest soln.

c  ## at later times want to use newest soln for spatial error flagging
c  ## at initial time want to use initial conditions (so retain symmetry for example)
         if (start_time+possk(lcheck) .ne. time) then  ! exact equality test here. counting on ieee arith.
             do 10 i = 1, mitot*mjtot*nvar
 10             alloc(locbig+i-1) = alloc(locnew+i-1)

             call bound(time,nvar,nghost,alloc(locbig),mitot,mjtot,mptr,
     1                  alloc(locaux),naux)
         else   ! boundary values already in locold
             locold = node(store2,mptr)
             do 11 i = 1, mitot*mjtot*nvar
 11             alloc(locbig+i-1) = alloc(locold+i-1)
         endif
c
c get user flags for refinement, which might be based on spatial gradient, 
c for example.  Use old values of soln at time t  before
c integration to get accurate boundary gradients
c
      if (flag_gradient) then
! need at least as big as nghost to fit ghost cells. if ibuff is bigger make
! the flagged array bigger so can buffer in place
         mbuff = max(nghost,ibuff+1)  
         mibuff = nx + 2*mbuff  !NOTE THIS NEW DIMENSIONING 
c                               !TO ALLOW ROOM FOR BUFFERING IN PLACE
         mjbuff = ny + 2*mbuff
         locamrflags = igetsp(mibuff*mjbuff)
         node(storeflags,mptr) = locamrflags

            do 20 i = 1, mibuff*mjbuff
 20         alloc(locamrflags+i-1) = goodpt

c        # call user-supplied routine to flag any points where 
c        # refinement is desired based on user's criterion.  
c        # Default version compares spatial gradient to tolsp.

          call flag2refine2(nx,ny,nghost,mbuff,nvar,naux,xleft,ybot,
     &        dx,dy,time,lcheck,tolsp,alloc(locbig),
     &        alloc(locaux),alloc(locamrflags),goodpt,badpt)
      endif

      end do
c      mptr = node(levelptr,mptr)
c      if (mptr .ne. 0) go to 5
c
      return
      end

