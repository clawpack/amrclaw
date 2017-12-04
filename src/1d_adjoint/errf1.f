c
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,rctcrse,mptr,mi2tot,
     2                 mitot,rctflg,mibuff,auxfine,
     2                 naux,nx,mbc)
      use amr_module
      use innerprod_module, only : calculate_innerproduct
      use adjoint_module
      implicit double precision (a-h,o-z)

 
      dimension  rctfine(nvar,mitot)
      dimension  auxfine(naux,mitot)
      dimension  aux_temp(naux,mitot)
      dimension  rctcrse(nvar,mi2tot)
      dimension  rctflg(mibuff)
      dimension  est(nvar,mitot)
      double precision levtol
      logical :: mask_selecta(totnum_adjoints)
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
      mask_selecta = .false.

c     Calculating correct tol for this level
      levtol = tol/numcells(1)
      do levcur = 1,levm-1
          levtol = levtol/intratx(levcur)
      enddo
 
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
      do  i  = nghost+1, mi2tot-nghost
          rflag = goodpt
          xofi  = xleft + (dble(ifine) - .5d0)*hx

c         calculate error in each term of coarse grid
          do k = 1,nvar
              term1 = rctfine(k,ifine)
              term2 = rctfine(k,ifine+1)
c             # divide by (aval*order) for relative error
              aval  = (term1+term2)/2.d0
              est(k,ifine)   =  dabs((aval-rctcrse(k,i))/ order)

c             retaining directionality of the wave
              est(k,ifine) = sign(est(k,ifine),rctcrse(k,i))

              est(k,ifine+1) = est(k,ifine)
          enddo

          ifine = ifine + 2
      enddo

      ! Loop over adjoint snapshots
      do k=1,totnum_adjoints

          ! Consider only snapshots that are within the desired time range
          if ((time+adjoints(k)%time) >= trange_start .and.
     .         (time+adjoints(k)%time) <= trange_final) then
                mask_selecta(k) = .true.
          endif
      enddo

      do k=1,totnum_adjoints-1
          if((.not. mask_selecta(k)) .and.
     .           (mask_selecta(k+1))) then
              mask_selecta(k) = .true.
              exit
          endif
      enddo

      do k=totnum_adjoints,2,-1
          if((.not. mask_selecta(k)) .and.
     .              (mask_selecta(k-1))) then
              mask_selecta(k) = .true.
              exit
          endif
      enddo

      ! Loop over adjoint snapshots
      do k=1,totnum_adjoints

          ! Consider only snapshots that are within the desired time range
          if (mask_selecta(k)) then

c             set innerproduct for fine grid
              aux_temp(innerprod_index,:) =
     .           calculate_innerproduct(time,est,k,nx,
     .           xleft,hx,nvar,mbc)

              ! Save max inner product
              do i=1,nx
                  auxfine(innerprod_index,i) =
     .              max(auxfine(innerprod_index,i),
     .              aux_temp(innerprod_index,i))
              enddo

          endif
      enddo

c     flag locations that need refining
      do i=1,nx

c          if (auxfine(innerprod_index,i) .gt. errmax)
c     .        errmax = auxfine(innerprod_index,i)
c          err2 = err2 + auxfine(innerprod_index,i)*
c     .                   auxfine(innerprod_index,i)

          if (auxfine(innerprod_index,i) .ge. levtol) then
               rflag  = badpt
          endif
          rctflg(i) = rflag
      enddo

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
