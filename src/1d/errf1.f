c
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,rctcrse,mptr,mi2tot,
     2                 mitot,rctflg,mibuff)
      use amr_module
      implicit double precision (a-h,o-z)

 
      dimension  rctfine(nvar,mitot)
      dimension  rctcrse(nvar,mi2tot)
      dimension  rctflg(mibuff)
c
c
c ::::::::::::::::::::::::::::: ERRF1 ::::::::::::::::::::::::::::::::
c
c  Richardson error estimator:  Used when flag_richardson is .true.
c  Compare error estimates in rctfine, rctcrse, 
c  A point is flagged if the error estimate is greater than tol
c  later we check if its in a region where its allowed to be flagged
c  or alternatively required.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

c
      time  = rnode(timemult, mptr)
      xleft = rnode(cornxlo,mptr)
      levm  = node(nestlevel, mptr)
      hx    = hxposs(levm)
      dt    = possk(levm)
      numsp = 0
 
      errmax = 0.0d0
      err2   = 0.0d0
c     order  = dt*dble(2**(iorder+1) - 2)
      order  = dble(2**(iorder+1) - 2)
c
      if (.not. (edebug)) go to 20
         write(outunit,107) mptr
 107     format(//,' coarsened grid values for grid ',i4)
         write(outunit,101) (rctcrse(1,i),
     .                          i = nghost+1, mi2tot-nghost)
         write(outunit,108) mptr
 108     format(//, ' fine grid values for grid ',i4)
         write(outunit,101) (rctfine(1,i),i=nghost+1,mitot-nghost)
101      format(' ',40e11.3)
c
c zero out the exterior locations so they don't affect err.est.
c
 20   continue
      ifine = nghost+1
c
      do 30  i  = nghost+1, mi2tot-nghost
        rflag = UNSET
c Only check errors if flag hasn't been set yet.
c If flag == DONTFLAG then refinement is forbidden by a region,
c if flag == DOFLAG checking is not needed
        if(rctflg(ifine) == UNSET
     .         .or. rctflg(ifine+1) == UNSET) then
          xofi  = xleft + (dble(ifine) - .5d0)*hx
          term1 = rctfine(1,ifine)
          term2 = rctfine(1,ifine+1)
c         # divide by (aval*order) for relative error
          aval  = (term1+term2)/2.d0
          est   =  dabs((aval-rctcrse(1,i))/ order)
          if (est .gt. errmax) errmax = est
          err2 = err2 + est*est
c         write(outunit,102) i,j,est,rctcrse(1,i,j)
 102      format(' i,est ',2i5,2e15.7)
c          write(outunit,104) term1,term2,term3,term4
 104      format('   ',4e15.7)
c         rctcrse(2,i,j) = est
c
          if (est .ge. tol) then
             rflag  = DOFLAG
          endif
        endif
      rctcrse(1,i) = rflag
      ifine = ifine + 2
 30   continue
c
c  print out intermediate flagged rctcrse (for debugging)
c
      if (eprint) then
         err2 = dsqrt(err2/dble((mi2tot-2*nghost)))
         write(outunit,103) mptr, levm, time,errmax, err2
 103     format(' grid ',i4,' level ',i4,' time ',e12.5,
     .          ' max. error = ',e15.7,' err2 = ',e15.7)
         if (edebug) then
           write(outunit,*) ' flagged points on coarsened grid ',
     .                      '(no ghost cells) for grid ',mptr
           write(outunit,106) (nint(rctcrse(1,i)),
     .                            i=nghost+1,mi2tot-nghost)
106           format(1h ,80i1)
         endif
      endif
c
      ifine   = nghost+1
      do 60 i = nghost+1, mi2tot-nghost
         if (rctcrse(1,i) .eq. DOFLAG) then
c           ## never set rctflg to DONTFLAG,
c           ## since flagregions or flag2refine may
c           ## have previously set it to DOFLAG
c           ## can only add DOFLAG pts in this routine
            rctflg(ifine)    = DOFLAG
            rctflg(ifine+1)  = DOFLAG
          endif
        ifine   = ifine + 2
 60     continue
c

      if (eprint) then
         write(outunit,118)
 118     format(' on fine grid (no ghost cells) flagged points are')
         if (edebug) then
           write(outunit,106)
     &      (nint(rctflg(i)),i=nghost+1,mitot-nghost)
        endif
      endif

      return
      end
