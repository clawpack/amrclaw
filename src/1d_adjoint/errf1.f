c
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,rctcrse,mptr,mi2tot,
     2                 mitot,rctflg,mibuff,auxfine,
     2                 naux,nx,mbc,jg)
      use amr_module
      use innerprod_module, only : calculate_innerproduct
      use adjoint_module, only: innerprod_index,
     .        totnum_adjoints, adjoints, trange_start, trange_final,
     .        levtol, eptr, errors
      implicit double precision (a-h,o-z)

 
      dimension  rctfine(nvar,mitot)
      dimension  auxfine(naux,mitot)
      dimension  aux_temp(naux,mitot)
      dimension  rctcrse(nvar,mi2tot)
      dimension  rctflg(mibuff)
      dimension  est(nvar,mitot)
      logical mask_selecta(totnum_adjoints), adjoints_found
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

      est(:,:) = 0.0d0
      auxfine(innerprod_index,:) = 0.0d0
      mask_selecta = .false.
      adjoints_found = .false.

c     Calculating correct tol for this level
c     --------------------
c     Total error allowed in this time step
      tol_exact = tol*dt/tfinal
c     Error allowed at this level
      tol_exact = tol_exact/(2**levm)
c     Error allowed per cell at this level
      tol_exact = tol_exact/(numcells(levm)*hx)

      if (t0+possk(levm) .eq. time) levtol(levm) = tol_exact
 
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
c         calculate error in each term of coarse grid
          do 50 k = 1,nvar
              term1 = rctfine(k,ifine)
              term2 = rctfine(k,ifine+1)
c             # divide by (aval*order) for relative error
              aval  = (term1+term2)/2.d0
              est(k,ifine)   =  dabs((aval-rctcrse(k,i))/ order)

c             retaining directionality of the wave
              est(k,ifine) = sign(est(k,ifine),rctcrse(k,i))
              est(k,ifine+1) = est(k,ifine)

 50       continue
      ifine = ifine + 2
 30   continue

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
              aux_temp(innerprod_index,:) =
     .              calculate_innerproduct(time,est,k,nx,
     .              xleft,hx,nvar,mbc)

              do 22  i  = 1, nx
                 auxfine(innerprod_index,i) =
     .              max(auxfine(innerprod_index,i),
     .              aux_temp(innerprod_index,i))
                 errors(eptr(jg)+i) = auxfine(innerprod_index,i)

                 if (auxfine(innerprod_index,i)
     .                               .ge. levtol(levm)) then
c                  ## never set rctflg to good, since flag2refine may
c                  ## have previously set it to bad
c                  ## can only add bad pts in this routine
                   rctflg(i)    = badpt
                 endif

 22           continue
          endif
 12   continue
c
c

      if (eprint) then
         write(outunit,118)
 118     format(' on fine grid (no ghost cells) flagged points are')
         if (edebug) then
           write(outunit,*)
     &      (nint(rctflg(i)),i=nghost+1,mitot-nghost)
        endif
      endif

      return
      end
