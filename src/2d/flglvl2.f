c
c -----------------------------------------------------------
c
      subroutine flglvl2(nvar,naux,lcheck,nxypts,index,lbase,npts,
     &                   start_time)
c
      use amr_module
      implicit double precision (a-h,o-z)
c
c :::::::::::::::::::: FLGLVL :::::::::::::::::::::::::::::::::
c
c flglvl = controls the error estimation/flagging bad pts. for
c          an entire level of grids.  returns pointer into alloc
c          where the (x,y) coordinations of the flagged pts. are.
c input parameters:
c           lcheck = level to be flagged
c output parameters:
c           nxypts = no. of flagged pts. total
c           index  = starting index in alloc of the flagged pts.
c                    (which occupy 2*nxypts locations).
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
      nxypts = 0
      numbad = 0

c     always call spest to set up stuff (initialize iflags, fill locbig)
      call spest2(nvar,naux,lcheck,start_time)
      if (flag_richardson) call errest(nvar,naux,lcheck)
c     if (tol .gt. 0.) call errest(nvar,naux,lcheck)

      call bufnst2(nvar,naux,numbad,lcheck,lbase) !NOW INCLUDES CALL TO DOMGRID INTERNALLY

      nxypts = nxypts + numbad
c
c  colate flagged pts into flagged points array
c  new version needs to check for proper nesting at this point
c  also needs to sort,  so can remove duplicates.
c
      if (nxypts .gt. 0) then  
c        build domain flags for each grid at level lcheck, instead of
c        previous approach using domain flags over entire  domain
c         call domgrid(lbase,lcheck)   ! will need since there are flagged pts  NOW IN BUFNST2
c
c in new version, there are bad cells but nxypts isnt true count any longer
c since there are duplicates, and proper nesting not yet checked
           index = igetsp(2*nxypts)
           call colate2(alloc(index),nxypts,lcheck,npts,lbase)
      else 
         npts = 0  !npts is number of unique flagged points after removing duplicates
         call freeFlags(lcheck)   ! otherwise storage freed in colate2. perhaps always do it here
      endif

      return
      end
c
c ---------------------------------------------------------------------------------
c
       subroutine freeFlags(lcheck)

       use amr_module
       implicit double precision (a-h, o-z)

       mptr = lstart(lcheck)
 10          continue
             locamrflags = node(storeflags,mptr)
             locdomflags = node(domflags_base,mptr)
             locdom2 = node(domflags2,mptr)

             nx     = node(ndihi,mptr) - node(ndilo,mptr) + 1
             ny     = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
             mbuff = max(nghost,ibuff+1)
             mibuff = nx + 2*mbuff
             mjbuff = ny + 2*mbuff

             ibytesPerDP = 8
             nwords = (mibuff*mjbuff)/ibytesPerDP+1
             call reclam(locdomflags, nwords)
             call reclam(locdom2, nwords)
             call reclam(locamrflags,mibuff*mjbuff)

        mptr = node(levelptr, mptr)
        if (mptr .ne. 0) go to 10

       return
       end
