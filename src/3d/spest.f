c
c -------------------------------------------------------------
c
      subroutine spest (nvar,naux,lcheck,dom1flags,isize,jsize,ksize,t0)
c
      use amr_module
      implicit double precision (a-h,o-z)

      integer*1  dom1flags (0:isize+1,0:jsize+1,0:ksize+1)

 
c :::::::::::::::::::::::::: SPEST :::::::::::::::::::::::::::::::::::
c for all grids at level lcheck:
c   Call user-supplied routine flag2refine to flag any points where
c   refinement is desired based on user's criterion.  
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c   initialize iflags here, so can put user flags right into it.
c
       do 4 i = 1, isize
       do 4 j = 1, jsize
       do 4 k = 1, ksize
 4        dom1flags(i,j,k) = 0
c
       mptr = lstart(lcheck)
 5     continue
          nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          nz     = node(ndkhi,mptr) - node(ndklo,mptr) + 1
          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          mktot  = nz + 2*nghost
          locnew = node(store1,mptr)
          locaux = node(storeaux,mptr)
          time   = rnode(timemult,mptr)
          dx     = hxposs(lcheck)
          dy     = hyposs(lcheck)
          dz     = hzposs(lcheck)
          xleft  = rnode(cornxlo,mptr)
          ybot   = rnode(cornylo,mptr)
          zbot   = rnode(cornzlo,mptr)
          xlow   = xleft - nghost*dx
          ylow   = ybot - nghost*dy
          zlow   = zbot - nghost*dz

c
          locbig = igetsp(mitot*mjtot*mktot*nvar)
          node(tempptr,mptr) = locbig
c         # straight copy into scratch array so don't mess up latest soln.
c  ## at later times want to use newest soln for spatial error flagging
c  ## at initial time want to use initial conditions (so retain symmetry for example)
!         if (t0+possk(lcheck) .ne. time) then  ! exact equality test here. counting on ieee arith.
            do 10 i = 1, mitot*mjtot*mktot*nvar
 10            alloc(locbig+i-1) = alloc(locnew+i-1)

            call bound(time,nvar,nghost,alloc(locbig),mitot,mjtot,mktot,
     1               mptr,alloc(locaux),naux)
!         else   ! boundary values already in locold
!             locold = node(store2,mptr)
!             do 11 i = 1, mitot*mjtot*nvar
! 11             alloc(locbig+i-1) = alloc(locold+i-1)
!         endif
c
c get user flags for refinement, which might be based on spatial gradient, 
c for example.  Use old values of soln at time t  before
c integration to get accurate boundary gradients
c
          if (flag_gradient) then
             locamrflags = igetsp(mitot*mjtot*mktot)
             do 20 i = 1, mitot*mjtot*mktot
 20             alloc(locamrflags+i-1) = goodpt
c
c        # call user-supplied routine to flag any points where 
c        # refinement is desired based on user's criterion.  
c        # Default version compares spatial gradient to tolsp.

         
         call flag2refine(nx,ny,nz,nghost,nvar,naux,xlow,ylow,zlow,
     &              dx,dy,dz,time,lcheck,tolsp,alloc(locbig),
     &              alloc(locaux),alloc(locamrflags), goodpt, badpt )
c

c        Put flags in iflags array now, so can reclaim space.
c        Note change of dimension of amrflags array:
c        3rd dim = 1 here, elsewhere flag array has idim3 = nvar
c
         idim3 = 1   ! 3rd dim = 1 here, elsewhere is nvar
         call setflags (dom1flags,isize,jsize,ksize,
     1                  alloc(locamrflags),idim3,mitot,mjtot,mktot,mptr)
         call reclam(locamrflags, mitot*mjtot*mktot)
      endif

c   previously used to reclam locbig space here. now save to reuse in errest
c   reclam locbig space afterwards.
c     call reclam(locbig,mitot*mjtot*nvar)

      mptr = node(levelptr,mptr)
      if (mptr .ne. 0) go to 5
c
      return
      end
