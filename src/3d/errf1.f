c
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,
     1                 rctcrse,mptr,mi2tot,mj2tot,mk2tot,
     2                 mitot,mjtot,mktot,rctflg)

      use amr_module
      implicit double precision (a-h,o-z)

 
      dimension  rctfine(nvar, mitot ,mjtot ,mktot )
      dimension  rctcrse(nvar, mi2tot,mj2tot,mk2tot)
      dimension  rctflg( nvar, mitot ,mjtot ,mktot )
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
            write(outunit,101) (rctcrse(1,i,j,k),
     .                          i = nghost+1, mi2tot-nghost)
10       continue
         write(outunit,108) mptr
 108     format(//, ' fine grid values for grid ',i4)
         do 15 kk = nghost+1, mktot-nghost
            k = mktot + 1 - kk
         do 15 jj = nghost+1, mjtot-nghost
            j = mjtot + 1 - jj
            write(outunit,101) (rctfine(1,i,j,k),
     .                          i = nghost+1, mitot-nghost)
15       continue
101      format(' ',13f6.3)
c
c zero out the exterior locations so they dont affect err.est.
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
          term1 = rctfine(1,ifine  ,jfine  ,kfine  )
          term2 = rctfine(1,ifine+1,jfine  ,kfine  )
          term3 = rctfine(1,ifine+1,jfine+1,kfine  )
          term4 = rctfine(1,ifine  ,jfine+1,kfine  )
          term5 = rctfine(1,ifine  ,jfine  ,kfine+1)
          term6 = rctfine(1,ifine+1,jfine  ,kfine+1)
          term7 = rctfine(1,ifine+1,jfine+1,kfine+1)
          term8 = rctfine(1,ifine  ,jfine+1,kfine+1)
c         # divide by (aval*order) for relative error
          aval  = (term1+term2+term3+term4
     &            +term5+term6+term7+term8)/8.d0
          est   =  dabs((aval-rctcrse(1,i,j,k))/ order)
          if (est .gt. errmax) errmax = est
            err2 = err2 + est*est
c         write(outunit,102) i,j,est
 102      format(' i,j,k,est ',3i5,e12.5)
c         rctcrse(2,i,j,k) = est
c
          if (est.ge.tol .and. allowflag(xofi,yofj,zofk,time,levm)) then
             rflag  = badpt
          endif 
      rctcrse(1,i,j,k) = rflag
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
 40      rctflg(1,i,j,k) = goodpt
c
c  print out intermediate flagged rctcrse (for debugging)
c
      if (eprint) then
         err2 = dsqrt(err2/dble((mi2tot-2*nghost)*(mj2tot-2*nghost)
     &                                           *(mk2tot-2*nghost)))
         write(outunit,103) mptr, levm, errmax, err2
 103     format(' grid ',i4,' level ',i4,
     .          ' max. error = ',e15.7,' err2 = ',e15.7)
         if (edebug) then
           write(outunit,*) ' flagged points on coarsened grid ',
     .                  'for grid ',mptr
           do 45 kk = nghost+1, mk2tot-nghost
              k = mk2tot + 1 - kk
           do 45 jj = nghost+1, mj2tot-nghost
              j = mj2tot + 1 - jj
              write(outunit,106) (nint(rctcrse(1,i,j,k)),
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
          if (rctcrse(1,i,j,k) .eq. goodpt) go to 55
            rctflg(1,ifine  ,jfine  ,kfine  ) = badpt
            rctflg(1,ifine+1,jfine  ,kfine  ) = badpt
            rctflg(1,ifine  ,jfine+1,kfine  ) = badpt
            rctflg(1,ifine+1,jfine+1,kfine  ) = badpt
            rctflg(1,ifine  ,jfine  ,kfine+1) = badpt
            rctflg(1,ifine+1,jfine  ,kfine+1) = badpt
            rctflg(1,ifine  ,jfine+1,kfine+1) = badpt
            rctflg(1,ifine+1,jfine+1,kfine+1) = badpt
 55   ifine   = ifine + 2
 60   continue
      jfine   = jfine + 2
 70   continue
      kfine   = kfine + 2
 74   continue
c
      return
      end
