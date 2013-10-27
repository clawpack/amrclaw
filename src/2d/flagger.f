c
c -----------------------------------------------------------
c
      subroutine flagger(nvar,naux,lcheck,start_time)

      use amr_module
      implicit double precision (a-h,o-z)

      integer omp_get_thread_num, omp_get_max_threads
      integer mythread/0/, maxthreads/1/
      integer listgrids(numgrids(lcheck))

c ::::::::::::::::::::: FLAGGER :::::::::::::::::::::::::
c
c flagger = set up for and call two routines that flag using
c    (a) spatial gradients, or other user-specified criteria
c    (b) richardson error estimates
c
c the two approaches share an array with boundary ghost values 
c
c ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      call prepgrids(listgrids,numgrids(lcheck),lcheck)
         mbuff = max(nghost,ibuff+1)  
c before parallel loop give grids the extra storage they need for error estimation
         do  jg = 1, numgrids(lcheck)
            mptr = listgrids(jg)
            nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
            ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
            mitot  = nx + 2*nghost
            mjtot  = ny + 2*nghost
            locbig = igetsp(mitot*mjtot*nvar)
            node(tempptr,mptr) = locbig
            mibuff = nx + 2*mbuff       ! NOTE THIS NEW DIMENSIONING 
            mjbuff = ny + 2*mbuff       ! TO ALLOW ROOM FOR BUFFERING IN PLACE
            locamrflags = igetsp(mibuff*mjbuff)  
            node(storeflags,mptr) = locamrflags
          end do

!$OMP PARALLEL DO PRIVATE(jg,mptr,nx,ny,mitot,mjtot,locnew,locaux),
!$OMP&            PRIVATE(time,dx,dy,xleft,ybot,xlow,ylow,locbig),
!$OMP&            PRIVATE(locold,mbuff,mibuff,mjbuff,locamrflags,i),
!$OMP&            SHARED(numgrids,listgrids,lcheck,nghost,nvar,naux),
!$OMP&            SHARED(tolsp,alloc,node,rnode,hxposs,hyposs,ibuff),
!$OMP&            SHARED(start_time,possk,flag_gradient,flag_richardson)
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
          locbig =  node(tempptr,mptr)
c         # straight copy into scratch array so don't mess up latest soln.

c         ## at later times want to use newest soln for spatial error flagging
c         ## at initial time want to use initial conditions (so retain symmetry for example)
         if (start_time+possk(lcheck) .ne. time) then !exact equality test-relying on ieee arith repeatability
c             do 10 i = 1, mitot*mjtot*nvar
c 10             alloc(locbig+i-1) = alloc(locnew+i-1)

c             call bound(time,nvar,nghost,alloc(locbig),mitot,mjtot,mptr,
c     1                  alloc(locaux),naux)


c            do in other order in case user messes up locbig in flag2refine, already have
c            them in locnew
             call bound(time,nvar,nghost,alloc(locnew),mitot,mjtot,mptr,
     1                  alloc(locaux),naux)
             do 10 i = 1, mitot*mjtot*nvar
 10             alloc(locbig+i-1) = alloc(locnew+i-1)

         else   ! boundary values already in locold
             locold = node(store2,mptr)
             do 11 i = 1, mitot*mjtot*nvar
 11             alloc(locbig+i-1) = alloc(locold+i-1)
         endif

!        # need at least as big as nghost to fit ghost cells. if ibuff is bigger make
!        # the flagged array bigger so can buffer in place
         mbuff = max(nghost,ibuff+1)  
         mibuff = nx + 2*mbuff       ! NOTE THIS NEW DIMENSIONING 
         mjbuff = ny + 2*mbuff       ! TO ALLOW ROOM FOR BUFFERING IN PLACE

!              ##  locamrflags used for flag storage. flag2refine flags directly into it.
!              ## richardson flags added to it. Then colate finished the job
                locamrflags = node(storeflags,mptr)
                do 20 i = 1, mibuff*mjbuff  ! initialize
 20                alloc(locamrflags+i-1) = goodpt

         if (flag_gradient) then

c     # call user-supplied routine to flag any points where 
c     # refinement is desired based on user's criterion.  
c     # Default version compares spatial gradient to tolsp.

            call flag2refine2(nx,ny,nghost,mbuff,nvar,naux,
     &                        xleft,ybot,dx,dy,time,lcheck,
     &                        tolsp,alloc(locbig),
     &                        alloc(locaux),alloc(locamrflags),
     &                        goodpt,badpt)
             endif     
c     
         if (flag_richardson) call errest(nvar,naux,lcheck,mptr,
     .                                    nx,ny)

       end do
! $OMP END PARALLEL DO

       return
       end
