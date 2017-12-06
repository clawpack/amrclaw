c
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,rctcrse,mptr,mi2tot,mj2tot,
     2                 mitot,mjtot,rctflg,mibuff,mjbuff,auxfine,
     2                 naux,nx,ny,mbc)
      use amr_module
      use innerprod_module, only : calculate_innerproduct
      use adjoint_module, only: innerprod_index,
     .        totnum_adjoints, adjoints, trange_start, trange_final
      implicit none

 
      real(kind=8),intent(in)::rctfine(nvar,1-mbc:nx+mbc,1-mbc:ny+mbc)
      real(kind=8),intent(in)::rctcrse(nvar,mi2tot,mj2tot)
      integer, intent(in) :: nx,ny,mbc,nvar,naux
      integer, intent(in) :: mibuff, mjbuff, mptr
      integer, intent(in) :: mitot,mjtot,mi2tot,mj2tot

      real(kind=8) :: rctflg(mibuff,mjbuff)
      real(kind=8) :: est(nvar,1-mbc:nx+mbc,1-mbc:ny+mbc)
      logical :: mask_selecta(totnum_adjoints)

      real(kind=8) :: auxfine(naux,1-mbc:nx+mbc,1-mbc:ny+mbc)
      real(kind=8) :: aux_temp(1:nx,1:ny)
      real(kind=8) :: term1,term2,term3,term4,aval

      real(kind=8) :: xleft,ybot,time,hx,hy,dt,order
      integer :: levm,i,j,ifine,jfine,jj,k

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
      auxfine(innerprod_index,:,:) = 0.0d0
      mask_selecta = .false.

c     order  = dt*dble(2**(iorder+1) - 2)
      order  = dble(2**(iorder+1) - 2)
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
      ifine = nghost+1
c
      do 30  i  = nghost+1, mi2tot-nghost
c         calculate error in each term of coarse grid
          do 50 k = 1,nvar
              term1 = rctfine(k,ifine,jfine)
              term2 = rctfine(k,ifine+1,jfine)
              term3 = rctfine(k,ifine+1,jfine+1)
              term4 = rctfine(k,ifine,jfine+1)
c             # divide by (aval*order) for relative error
              aval  = (term1+term2+term3+term4)/4.d0
              est(k,ifine,jfine)   =  dabs((aval-rctcrse(k,i,j))/ order)

c             retaining directionality of the wave
              est(k,ifine,jfine) =
     .            sign(est(k,ifine,jfine),rctcrse(k,i,j))
              est(k,ifine + 1,jfine) = est(k,ifine,jfine)
              est(k,ifine + 1,jfine + 1) = est(k,ifine,jfine)
              est(k,ifine,jfine + 1) = est(k,ifine,jfine)
 50       continue
      ifine = ifine + 2
 30   continue
      jfine = jfine + 2
 35   continue

c     Loop over adjoint snapshots
      do k=1,totnum_adjoints
          if ((time+adjoints(k)%time) >= trange_start .and.
     .        (time+adjoints(k)%time) <= trange_final) then
              mask_selecta(k) = .true.
          endif
      enddo

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
              aux_temp(:,:) =
     .            calculate_innerproduct(time,est,k,nx,ny,
     .            xleft,ybot,hx,hy,nvar,mbc)

              do 22  i  = 1, nx
              do 23  j  = 1, ny
                  auxfine(innerprod_index,i,j) =
     .                max(auxfine(innerprod_index,i,j),
     .                    aux_temp(i,j))

                  if (auxfine(innerprod_index,i,j) .ge. tol) then
c                     ## never set rctflg to good, since flag2refine may
c                     ## have previously set it to bad
c                     ## can only add bad pts in this routine
                      rctflg(i,j)    = badpt
                  endif

 23           continue
 22           continue
          endif
 12   continue
c

      if (eprint) then
         write(outunit,118)
 118     format(' on fine grid (no ghost cells) flagged points are')
         if (edebug) then
          do 56 jj = nghost+1, mjtot-nghost
           j = mjtot + 1 - jj
106           format(1h ,80i1)
           write(outunit,106)
     &      (nint(rctflg(i,j)),i=nghost+1,mitot-nghost)
 56       continue
        endif
      endif

      return
      end
