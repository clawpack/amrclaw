c
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,rctcrse,mptr,mi2tot,mj2tot,
     2                 mitot,mjtot,rctflg,mibuff,mjbuff,auxfine,
     2                 naux,nx,ny,mbc,jg)
      use amr_module
      use innerprod_module, only : calculate_innerproduct
      use adjoint_module, only: innerprod_index,
     .        totnum_adjoints, adjoints, trange_start, trange_final,
     .        levtol, eptr, errors
      implicit double precision (a-h,o-z)

 
      dimension  rctfine(nvar,mitot,mjtot)
      dimension  rctcrse(nvar,mi2tot,mj2tot)
      dimension  rctflg(mibuff,mjbuff)

      dimension  aux_crse(mi2tot,mj2tot)
      dimension  aux_temp(mi2tot,mj2tot)
      dimension  err_crse(nvar,mi2tot,mj2tot)
      logical mask_selecta(totnum_adjoints), adjoints_found
      dimension  auxfine(naux,mitot,mjtot)
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
      ybot  = rnode(cornylo,mptr)
      hy    = hyposs(levm)
      dt    = possk(levm)
      numsp = 0
 
      errmax = 0.0d0
      err2   = 0.0d0

      mask_selecta = .false.
      adjoints_found = .false.
      aux_crse = 0.0d0
      aux_temp = 0.0d0

c     order  = dt*dble(2**(iorder+1) - 2)
      order  = dble(2**(iorder+1) - 2)
c
c     Calculating correct tol for this level
c     nxnest is the maximum number of refinement levels, from amr_module
c     --------------------
      tol_exact = tol*dt/(tfinal-t0)
      tol_exact = tol_exact/mxnest
      tol_exact = tol_exact/(numcells(levm)*hx*hy)

      if (t0+possk(levm) .eq. time) levtol(levm) = tol_exact
c
      if (.not. (edebug)) go to 20
         write(outunit,107) mptr
 107     format(//,' coarsened grid values for grid ',i4)
         do 10 jj = nghost+1, mj2tot-nghost
            j = mj2tot + 1 - jj
            write(outunit,101) (rctcrse(1,i,j),
     .                          i = nghost+1, mi2tot-nghost)
10       continue
         write(outunit,108) mptr
 108     format(//, ' fine grid values for grid ',i4)
         do 15 jj = nghost+1, mjtot-nghost
            j = mjtot + 1 - jj
            write(outunit,101) (rctfine(1,i,j),i=nghost+1,mitot-nghost)
15       continue
101      format(' ',40e11.3)
c
c zero out the exterior locations so they don't affect err.est.
c
 20   continue
      jfine = nghost+1
      do 35  j = nghost+1, mj2tot-nghost
      yofj  = ybot + (dble(jfine) - .5d0)*hy
      ifine = nghost+1
c
      do 30  i  = nghost+1, mi2tot-nghost
        do 40 k = 1, nvar
          xofi  = xleft + (dble(ifine) - .5d0)*hx
          term1 = rctfine(k,ifine,jfine)
          term2 = rctfine(k,ifine+1,jfine)
          term3 = rctfine(k,ifine+1,jfine+1)
          term4 = rctfine(k,ifine,jfine+1)
c         # divide by (aval*order) for relative error
          aval  = (term1+term2+term3+term4)/4.d0
          est   =  dabs((aval-rctcrse(k,i,j))/ order)
          if (est .gt. errmax) errmax = est
          err2 = err2 + est*est
c         write(outunit,102) i,j,est,rctcrse(1,i,j)
 102      format(' i,j,est ',2i5,2e15.7)
c          write(outunit,104) term1,term2,term3,term4
 104      format('   ',4e15.7)
c         rctcrse(2,i,j) = est
c
          err_crse(k,i,j) = est
c         retaining directionality of the wave
          err_crse(k,i,j) = sign(est,rctcrse(k,i,j))
 40     continue
      ifine = ifine + 2
 30   continue
      jfine = jfine + 2
 35   continue

c     Loop over adjoint snapshots
      do k=1,totnum_adjoints
          if ((time+adjoints(k)%time) >= trange_start .and.
     .        (time+adjoints(k)%time) <= trange_final) then
              mask_selecta(k) = .true.
              adjoints_found = .true.
          endif
      enddo

      if(.not. adjoints_found) then
          write(*,*) "Error: no adjoint snapshots ",
     .        "found in time range."
          write(*,*) "Consider increasing time rage of interest, ",
     .        "or adding more snapshots."
      endif

      do k=1,totnum_adjoints-1
          if((.not. mask_selecta(k)) .and.
     .        (mask_selecta(k+1))) then
              mask_selecta(k) = .true.
              exit
          endif
      enddo

      do k=totnum_adjoints,2,-1
          if((.not. mask_selecta(k)) .and.
     .        (mask_selecta(k-1))) then
              mask_selecta(k) = .true.
              exit
          endif
      enddo

      do 12 k = 1,totnum_adjoints
c         ! Consider only snapshots that are within the desired time range
          if (mask_selecta(k)) then

c             set innerproduct for fine grid
              aux_temp(
     .              nghost+1:mi2tot-nghost, nghost+1:mj2tot-nghost) =
     .              calculate_innerproduct(time,err_crse,k,nx/2,
     .              ny/2,xleft,ybot,hx*2,hy*2,nvar,mbc)

              do 22  i  = nghost+1, mi2tot-nghost
                do 32  j = nghost+1, mj2tot-nghost
                  aux_crse(i,j) = max(aux_crse(i,j),aux_temp(i,j))

                  i_val = i-nghost
                  j_val = j-nghost
                  errors(eptr(jg)+(i_val-1)*ny/2+j_val) = aux_crse(i,j)

                  rctcrse(1,i,j)  = goodpt
                  if (aux_crse(i,j) .ge. levtol(levm)) then
c                    ## never set rctflg to good, since flag2refine may
c                    ## have previously set it to bad
c                    ## can only add bad pts in this routine
                      rctcrse(1,i,j)  = badpt
                  endif

 32             continue
 22           continue
          endif
 12   continue
c
c  print out intermediate flagged rctcrse (for debugging)
c
      if (eprint) then
         err2 = dsqrt(err2/dble((mi2tot-2*nghost)*(mj2tot-2*nghost)))
         write(outunit,103) mptr, levm, time,errmax, err2
 103     format(' grid ',i4,' level ',i4,' time ',e12.5,
     .          ' max. error = ',e15.7,' err2 = ',e15.7)
         if (edebug) then
           write(outunit,*) ' flagged points on coarsened grid ',
     .                      '(no ghost cells) for grid ',mptr
           do 45 jj = nghost+1, mj2tot-nghost
              j = mj2tot + 1 - jj
              write(outunit,106) (nint(rctcrse(1,i,j)),
     .                            i=nghost+1,mi2tot-nghost)
106           format(1h ,80i1)
45         continue
         endif
      endif
c
      jfine   = nghost+1
      do 70 j = nghost+1, mj2tot-nghost
      ifine   = nghost+1
      do 60 i = nghost+1, mi2tot-nghost
         auxfine(innerprod_index,ifine,jfine) = aux_crse(i,j)
         auxfine(innerprod_index,ifine+1,jfine) = aux_crse(i,j)
         auxfine(innerprod_index,ifine,jfine+1) = aux_crse(i,j)
         auxfine(innerprod_index,ifine+1,jfine+1) = aux_crse(i,j)

         if (rctcrse(1,i,j) .eq. goodpt) go to 55
c           ## never set rctflg to good, since flag2refine may
c           ## have previously set it to bad
c           ## can only add bad pts in this routine
            rctflg(ifine,jfine)    = badpt
            rctflg(ifine+1,jfine)  = badpt
            rctflg(ifine,jfine+1)  = badpt
            rctflg(ifine+1,jfine+1)= badpt
 55       ifine   = ifine + 2
 60     continue
        jfine   = jfine + 2
 70   continue
c

      if (eprint) then
         write(outunit,118)
 118     format(' on fine grid (no ghost cells) flagged points are')
         if (edebug) then
          do 56 jj = nghost+1, mjtot-nghost
           j = mjtot + 1 - jj
           write(outunit,106)
     &      (nint(rctflg(i,j)),i=nghost+1,mitot-nghost)
 56       continue
        endif
      endif

      return
      end
