c
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,
     1                 rctcrse,mptr,mi2tot,mj2tot,mk2tot,
     2                 mitot,mjtot,mktot,rctflg)
      implicit double precision (a-h,o-z)

      include  "call.i"
 
      dimension  rctfine(mitot ,mjtot ,mktot ,nvar)
      dimension  rctcrse(mi2tot,mj2tot,mk2tot,nvar)
      dimension  rctflg( mitot ,mjtot ,mktot ,nvar)
      logical    allowflag
      external   allowflag
c
c :::::::::::::::::::: ERRF1 :::::::::::::::::::::::::::::::::::::
c
c  Richardson error estimator: used when tol>0 in user input.
c  Compare error estimates in rctfine, rctcrse.
c  A point is flagged if the error estimate is greater than tol
c  and is allowflag(x,y,t,level) = .true. at this point.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      time  = rnode(timemult, mptr)
      xleft = rnode(cornxlo,mptr)
      levm  = node(nestlevel, mptr)
      hx    = hxposs(levm)
      yfront= rnode(cornylo,mptr)
      hy    = hyposs(levm)
      zbot  = rnode(cornzlo,mptr)
      hz    = hzposs(levm)
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
         do 10 kk = nghost+1, mk2tot-nghost
            k = mk2tot + 1 - kk
         do 10 jj = nghost+1, mj2tot-nghost
            j = mj2tot + 1 - jj
            write(outunit,101) (rctcrse(i,j,k,1),
     .                          i = nghost+1, mi2tot-nghost)
10       continue
         write(outunit,108) mptr
 108     format(//, ' fine grid values for grid ',i4)
         do 15 kk = nghost+1, mktot-nghost
            k = mktot + 1 - kk
         do 15 jj = nghost+1, mjtot-nghost
            j = mjtot + 1 - jj
            write(outunit,101) (rctfine(i,j,k,1),
     .                          i = nghost+1, mitot-nghost)
15       continue
101      format(' ',13f6.3)
c
c zero out the exterior locations so they don't affect err.est.
c
 20   continue
      kfine = nghost+1
      do 39  k = nghost+1, mk2tot-nghost
      zofk  = zbot   + (dble(kfine) - .5d0)*hz
      jfine = nghost+1
      do 35  j = nghost+1, mj2tot-nghost
      yofj  = yfront + (dble(jfine) - .5d0)*hy
      ifine = nghost+1
c
      do 30  i  = nghost+1, mi2tot-nghost
          rflag = goodpt
          xofi  = xleft + (dble(ifine) - .5d0)*hx
          term1 = rctfine(ifine  ,jfine  ,kfine  ,1)
          term2 = rctfine(ifine+1,jfine  ,kfine  ,1)
          term3 = rctfine(ifine+1,jfine+1,kfine  ,1)
          term4 = rctfine(ifine  ,jfine+1,kfine  ,1)
          term5 = rctfine(ifine  ,jfine  ,kfine+1,1)
          term6 = rctfine(ifine+1,jfine  ,kfine+1,1)
          term7 = rctfine(ifine+1,jfine+1,kfine+1,1)
          term8 = rctfine(ifine  ,jfine+1,kfine+1,1)
c         # divide by (aval*order) for relative error
          aval  = (term1+term2+term3+term4
     &            +term5+term6+term7+term8)/8.d0
          est   =  dabs((aval-rctcrse(i,j,k,1))/ order)
          if (est .gt. errmax) errmax = est
	        err2 = err2 + est*est
c         write(outunit,102) i,j,est
 102      format(' i,j,est ',2i5,e12.5)
c         rctcrse(i,j,2) = est
c
          if (est.ge.tol .and. allowflag(xofi,yofj,zofk,time,levm)) then
             rflag  = badpt
          endif 
      rctcrse(i,j,k,1) = rflag
      ifine = ifine + 2
 30   continue
      jfine = jfine + 2
 35   continue
      kfine = kfine + 2
 39   continue
c
c  transfer flagged points on cell centered coarse grid
c  to cell centered fine grid. count flagged points.
c
c  initialize rctflg to 0.0 (no flags)  before flagging
c
      do 40 k = 1, mktot
      do 40 j = 1, mjtot
      do 40 i = 1, mitot
 40      rctflg(i,j,k,1) = goodpt
c
c  print out intermediate flagged rctcrse (for debugging)
c
      if (eprint) then
	 err2 = dsqrt(err2/dble((mi2tot-2*nghost)*(mj2tot-2*nghost)
     &                                             *(mk2tot-2*nghost)))
         write(outunit,103) mptr, levm, errmax, err2
 103     format(' grid ',i4,' level ',i4,
     .          ' max. error = ',e15.7,' err2 = ',e15.7)
	 if (edebug) then
	   write(outunit,*) ' flagged points on coarsened grid ',
     .             	    'for grid ',mptr
           do 45 kk = nghost+1, mk2tot-nghost
              k = mk2tot + 1 - kk
	   do 45 jj = nghost+1, mj2tot-nghost
	      j = mj2tot + 1 - jj
	      write(outunit,106) (nint(rctcrse(i,j,k,1)),
     .                            i=nghost+1,mi2tot-nghost)
106           format(1h ,80i1)
45         continue
         endif
      endif
c
      kfine   = nghost+1
      do 74 k = nghost+1, mk2tot-nghost
      jfine   = nghost+1
      do 70 j = nghost+1, mj2tot-nghost
      ifine   = nghost+1
      do 60 i = nghost+1, mi2tot-nghost
	 if (rctcrse(i,j,k,1) .eq. goodpt) go to 55
	    rctflg(ifine  ,jfine  ,kfine  ,1) = badpt
	    rctflg(ifine+1,jfine  ,kfine  ,1) = badpt
	    rctflg(ifine  ,jfine+1,kfine  ,1) = badpt
	    rctflg(ifine+1,jfine+1,kfine  ,1) = badpt
            rctflg(ifine  ,jfine  ,kfine+1,1) = badpt
            rctflg(ifine+1,jfine  ,kfine+1,1) = badpt
            rctflg(ifine  ,jfine+1,kfine+1,1) = badpt
            rctflg(ifine+1,jfine+1,kfine+1,1) = badpt
 55   ifine   = ifine + 2
 60   continue
      jfine   = jfine + 2
 70   continue
      kfine   = kfine + 2
 74   continue
c
c CHANGED **********************
c ** spatial errore stimated in sperr.f now using user routine flag2refine.f *****
c      if (edebug) then
c	       write(outunit,*)" spatial error for grid ",mptr
c         do 75 kkfine = nghost+1, mktot-nghost
c            kfine = mktot + 1 - kkfine
c         do 75 jjfine = nghost+1, mjtot-nghost
c	          jfine = mjtot + 1 - jjfine
c	          write(outunit,101)(sperr(ifine,jfine,kfine),
c     .		          	       ifine=nghost+1,mitot-nghost)
c 75      continue
c     endif
c
c      do 80 kfine = nghost+1, mktot-nghost
c      zofk  = zbot + (dble(kfine) - nghost - .5d0)*hz
c      do 80 jfine = nghost+1, mjtot-nghost
c      yofj  = yfront + (dble(jfine) - nghost - .5d0)*hy
c      do 80 ifine = nghost+1, mitot-nghost
c        xofi  = xleft + (dble(ifine) - nghost - .5d0)*hx
c        if (sperr(ifine,jfine,kfine) .gt. tolsp .and.
c     &      allowed(xofi,yofj,zofk,levm)) then
c	       rflag = rctflg(ifine,jfine,kfine,1)
c               if (rflag .ne. badpt) then
c	         rctflg(ifine,jfine,kfine,1) = badpt
c      	         numsp = numsp + 1
c               endif
c	  endif
c 80   continue

c      if (eprint) then
c	 write(outunit,118) numsp,mptr
c 118     format( i5,' more pts. flagged for spatial error on grid',i4,/)
c	if (edebug) then
c          do 56 kk = nghost+1, mktot-nghost
c          k = mktot + 1 - kk
c          write(6,*) "z plane = ",k
c
c	  do 56 jj = nghost+1, mjtot-nghost
c	   j = mjtot + 1 - jj
c	   write(outunit,106)(nint(rctflg(i,j,k,1)),
c     &                                    i=nghost+1,mitot-nghost)
c 56       continue
c	endif
c
c      endif
c END OF CHANGE *******************************

      return
      end
