! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::     Support module for the adjoint_module.
! :::::: This module contains subroutines that are specific 
! :::::: to AMRClaw, and is mirrored by a similar module that
! :::::: is specific to GeoClaw.
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module adjointsup_module

contains

! ========================================================================
!  Routine to calculate inner product
! ========================================================================

    subroutine calculate_innerproduct(q,k,mx_f,my_f,xlower_f, &
               ylower_f,dx_f,dy_f,meqn_f,mbc_f,maux_f,aux)

        use adjoint_module

        implicit none

        real(kind=8), intent(in) :: xlower_f,ylower_f,dx_f,dy_f
        integer :: k,mx_f,my_f,meqn_f,mbc_f,maux_f
        real(kind=8), intent(in) :: q(meqn_f,1-mbc_f:mx_f+mbc_f,1-mbc_f:my_f+mbc_f)

        integer :: mx_a, my_a, mptr_a, mbc_a
        integer :: i, j, i1, i2, j1, j2, level, loc, z
        real(kind=8) :: dx_a, xlower_a, xupper_a, xupper_f
        real(kind=8) :: dy_a, ylower_a, yupper_a, yupper_f
        real(kind=8) :: x1, x2, y1, y2

        real(kind=8) :: q_innerprod(mx_f,my_f)
        logical :: mask_forward(mx_f,my_f)
        real(kind=8) :: q_interp(meqn_f,mx_f,my_f)
        real(kind=8), intent(inout) :: aux(maux_f,1-mbc_f:mx_f+mbc_f,1-mbc_f:my_f+mbc_f)

        logical, allocatable :: mask_adjoint(:,:)

        xupper_f = xlower_f + mx_f*dx_f
        yupper_f = ylower_f + my_f*dy_f

        ! Loop over patches in adjoint solution
        do z = 1, adjoints(k)%ngrids
            mptr_a = adjoints(k)%gridpointer(z)
            level = adjoints(k)%gridlevel(mptr_a)

            ! Number of points in x and y
            mx_a = adjoints(k)%ncellsx(mptr_a)
            my_a = adjoints(k)%ncellsy(mptr_a)

            ! Finding extreem values for grid
            xlower_a = adjoints(k)%xlowvals(mptr_a)
            dx_a = adjoints(k)%hxposs(level)
            xupper_a = xlower_a + mx_a*dx_a
            ylower_a = adjoints(k)%ylowvals(mptr_a)
            dy_a = adjoints(k)%hyposs(level)
            yupper_a = ylower_a + my_a*dy_a

            loc = adjoints(k)%loc(mptr_a)

            ! Check if adjoint patch overlaps with forward patch
            x1 = max(xlower_f,xlower_a)
            x2 = min(xupper_f,xupper_a)
            y1 = max(ylower_f,ylower_a)
            y2 = min(yupper_f,yupper_a)

            if ((x1 > x2) .or. (y1 > y2)) then
                ! Skipping interpolation if grids don't overlap
                mask_forward = .false.
                continue
            else
                mbc_a = adjoints(k)%nghost
                allocate(mask_adjoint(1-mbc_a:mx_a+mbc_a, 1-mbc_a:my_a+mbc_a))

                ! Create a mask that is .true. only in part of patch intersecting forward patch:
                i1 = max(int((x1 - xlower_a + 0.5d0*dx_a) / dx_a), 0)
                i2 = min(int((x2 - xlower_a + 0.5d0*dx_a) / dx_a) + 1, mx_a+1)
                j1 = max(int((y1 - ylower_a + 0.5d0*dy_a) / dy_a), 0)
                j2 = min(int((y2 - ylower_a + 0.5d0*dy_a) / dy_a) + 1, my_a+1)

                forall (i=1-mbc_a:mx_a+mbc_a, j=1-mbc_a:my_a+mbc_a)
                    mask_adjoint(i,j) = ((i >= i1) .and. (i <= i2) .and. &
                                       (j >= j1) .and. (j <= j2))
                end forall

                ! Create a mask that is .true. only in part of forward patch intersecting patch:

                i1 = max(int((x1 - xlower_f + 0.5d0*dx_f) / dx_f)+1, 0)
                i2 = min(int((x2 - xlower_f + 0.5d0*dx_f) / dx_f), mx_f)
                j1 = max(int((y1 - ylower_f + 0.5d0*dy_f) / dy_f)+1, 0)
                j2 = min(int((y2 - ylower_f + 0.5d0*dy_f) / dy_f), my_f)

                forall (i=1:mx_f, j=1:my_f)
                    mask_forward(i,j) = ((i >= i1) .and. (i <= i2) .and. &
                                       (j >= j1) .and. (j <= j2))
                end forall

                ! Interpolate adjoint values to q_interp
                ! Note that values in q_interp will only be set properly where
                ! mask_adjoint == .true.
                call interp_adjoint( &
                        adjoints(k)%meqn, k, q_interp, xlower_a, ylower_a, &
                        dx_a, dy_a, mx_a, my_a, xlower_f, ylower_f, &
                        dx_f, dy_f, mx_f, my_f, &
                        mask_adjoint, mptr_a, mask_forward)

                q_innerprod = 0.d0

                ! For each overlapping point, calculate inner product
                ! When using richardson error, q contains the error 
                ! in each term.
                forall(i = 1:mx_f, j = 1:my_f, mask_forward(i,j))
                    q_innerprod(i,j) = abs(dot_product(q(:,i,j),q_interp(:,i,j)))
                end forall

                do i=1,mx_f
                    do j=1,my_f
                        if (q_innerprod(i,j) > aux(innerprod_index,i,j)) then
                            aux(innerprod_index,i,j) = q_innerprod(i,j)
                        endif
                    enddo
                enddo

                deallocate(mask_adjoint)
            endif
        enddo

    end subroutine calculate_innerproduct

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine interp_adjoint(nvar, r, q_interp, xlower_a, ylower_a, &
              dx_a, dy_a, mx_a, my_a, xlower_f, ylower_f, dx_f, &
              dy_f, mx_f, my_f, mask_adjoint, mptr_a, mask_forward)

        use adjoint_module, only : adjoints
        implicit none

        ! Function arguments
        integer, intent(in) :: nvar, r, mx_a, my_a
        logical, intent(in) :: mask_adjoint( &
                          1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost, &
                          1-adjoints(r)%nghost:my_a+adjoints(r)%nghost)
        real(kind=8), intent(in) :: xlower_a, xlower_f, ylower_a, ylower_f
        integer, intent(in) :: mx_f, my_f, mptr_a
        real(kind=8), intent(in) :: dx_f, dx_a, dy_f, dy_a

        integer :: z, k, iz, jk, mitot
        integer :: ivar, i, j, iadd, loc
        real(kind=8) :: q_interp(nvar,mx_f,my_f), denom
        real(kind=8) :: x, xhigh_a, y, yhigh_a
        real(kind=8) :: dxz, dyk, a, b, c
        logical :: mask_forward(mx_f,my_f)

        iadd(ivar,i,j)  = loc + ivar - 1 + adjoints(r)%meqn*((j-1)*mitot+i-1)

        q_interp = 0.0
        xhigh_a  = xlower_a + mx_a*dx_a
        yhigh_a = ylower_a + my_a*dx_a
        loc    = adjoints(r)%loc(mptr_a)
        mitot = adjoints(r)%ncellsx(mptr_a) + 2*adjoints(r)%nghost

        do z = 1,mx_f
            do k = 1,my_f
                if (mask_forward(z,k)) then
                    x = xlower_f + (z - 0.5d0)*dx_f
                    y = ylower_f + (k - 0.5d0)*dy_f

                    iz = int((x - xlower_a + 0.5d0*dx_a) / dx_a) + 1
                    dxz = x - (xlower_a + (iz-0.5d0)*dx_a)
                    jk = int((y - ylower_a + 0.5d0*dy_a) / dy_a) + 1
                    dyk = y - (ylower_a + (jk-0.5d0)*dy_a)

                    ! Interpolate only if this cell is overlapping with grid
                    if (mask_adjoint(iz,jk)) then
                        do ivar=1,nvar

                        a = (adjoints(r)%alloc(iadd(ivar,iz+1,jk)) &
                              - adjoints(r)%alloc(iadd(ivar,iz,jk))) / dx_a
                        b = (adjoints(r)%alloc(iadd(ivar,iz,jk+1)) &
                              - adjoints(r)%alloc(iadd(ivar,iz,jk))) / dy_a
                        c = (adjoints(r)%alloc(iadd(ivar,iz+1,jk+1)) &
                              + adjoints(r)%alloc(iadd(ivar,iz,jk)) &
                              - adjoints(r)%alloc(iadd(ivar,iz+1,jk)) &
                              - adjoints(r)%alloc(iadd(ivar,iz,jk+1))) / (dx_a*dy_a)

                        q_interp(ivar,z,k) = &
                              adjoints(r)%alloc(iadd(ivar,iz,jk)) &
                              + a*dxz + b*dyk + c*dxz*dyk

                        enddo
                    endif
                endif
            enddo
        enddo

    end subroutine interp_adjoint

! ========================================================================
!  Routine to compute error in forward solution for the adjoint method
!  Adapted from errf1.f
!  Differences: 
!   (i) Computes error for all of the terms in q
!   (ii) Needs the aux array passed in, because the inner product is saved 
!        to the aux array (note, this changes the call sequence, which is 
!        why this is not included into the AMRClaw errf1.f file. If plotting 
!        of the inner product is never wanted, this could be simplified 
!        and included into the AMRClaw errf1.f file.
! ========================================================================

    subroutine errf1a(rctfine,nvar,rctcrse,mptr,mi2tot,mj2tot, &
                      mitot,mjtot,rctflg,mibuff,mjbuff,auxfine, &
                      naux,auxcrse)
      use amr_module
      use adjoint_module, only: innerprod_index, &
             totnum_adjoints, adjoints, trange_start, trange_final, &
             levtol, eptr, errors, grid_num, select_snapshots
      implicit double precision (a-h,o-z)

      dimension  rctfine(nvar,mitot,mjtot)
      dimension  rctcrse(nvar,mi2tot,mj2tot)
      dimension  rctflg(mibuff,mjbuff)
      dimension  err_crse(nvar,mi2tot,mj2tot)
      dimension  bcrse(mitot,mjtot)
      dimension  auxfine(naux,mitot,mjtot)
      dimension  auxcrse(naux,mi2tot,mj2tot)

      logical mask_selecta(totnum_adjoints)
!
!
! :::::::::::::::::::::: Modified from ERRF1 :::::::::::::::::::::::::
!
!  Richardson error estimator:  Used when flag_richardson is .true.
!  Compare error estimates in rctfine, rctcrse,
!  A point is flagged if the inner product between the
!  error estimate and the adjoint solution is greater than tol
!  later we check if its in a region where its allowed to be flagged
!  or alternatively required.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!

      time  = rnode(timemult, mptr)
      xleft = rnode(cornxlo,mptr)
      levm  = node(nestlevel, mptr)
      hx    = hxposs(levm)
      ybot  = rnode(cornylo,mptr)
      hy    = hyposs(levm)
      dt    = possk(levm)

      jg = grid_num(mptr)
      nx = mitot - 2*nghost
      ny = mjtot - 2*nghost

      call select_snapshots(time,mask_selecta)

      errmax = 0.0d0
      err2   = 0.0d0
      auxfine(innerprod_index,:,:) = 0.0d0
      auxcrse(innerprod_index,:,:) = 0.0d0

      order  = dble(2**(iorder+1) - 2)
!
!     Calculating correct tol for this level
!     nxnest is the maximum number of refinement levels, from amr_module
!     --------------------
!     Total error allowed in this time step
      tol_exact = tol*dt/(tfinal-t0)
!     Error allowed at this level
      tol_exact = tol_exact/mxnest
!     Error allowed per cell at this level
      tol_exact = tol_exact/(numcells(levm)*hx*hy)

      if (t0+possk(levm) .eq. time) levtol(levm) = tol_exact

!
! Main loop: calculate errors in forward problem
!
      jfine = nghost+1
      do j = nghost+1, mj2tot-nghost
        yofj  = ybot + (dble(jfine) - .5d0)*hy
        ifine = nghost+1

        do i  = nghost+1, mi2tot-nghost

! Only check errors if flag hasn't been set yet.
! If flag == DONTFLAG then refinement is forbidden by a region,
! if flag == DOFLAG checking is not needed

! Note: here rctcrse is being used as a temporary flag
! the fine grid amrflags array is stored in rctflg, and will be
! updated based on rctcrse at the end of this routine
           if(rctflg(ifine,jfine) == UNSET &
              .or. rctflg(ifine+1,jfine) == UNSET &
              .or. rctflg(ifine,jfine+1) == UNSET &
              .or. rctflg(ifine+1,jfine+1) == UNSET) then

              xofi  = xleft + (dble(ifine) - .5d0)*hx

!             Calculate error for each value of q

              do k = 1, nvar
                term1 = rctfine(k,ifine,jfine)
                term2 = rctfine(k,ifine+1,jfine)
                term3 = rctfine(k,ifine+1,jfine+1)
                term4 = rctfine(k,ifine,jfine+1)
!               # divide by (aval*order) for relative error
                aval  = (term1+term2+term3+term4)/4.d0
                est   =  dabs((aval-rctcrse(k,i,j))/ order)
                if (est .gt. errmax) errmax = est
                  err2 = err2 + est*est

                  err_crse(k,i,j) = est
!                 retaining directionality of the wave
                  err_crse(k,i,j) = sign(est,rctcrse(k,i,j))
              enddo
            else
                err_crse(:,i,j) = 0.d0
            endif

          ifine = ifine + 2
        enddo

        jfine = jfine + 2
      enddo

!     Loop over adjoint snapshots

! Note: here the coarse version of the inner product values in auxcrse are set.
! The fine grid inner product values will be stored in 
! auxfine(inneprod_index,:,:) which will be updated based 
! on rctcrse at the end of this routine. 

      do k = 1,totnum_adjoints
!         Consider only snapshots that are within the desired time range
          if (mask_selecta(k)) then
!             set innerproduct
              call calculate_innerproduct(err_crse,k,nx/2, &
                    ny/2,xleft,ybot,hx*2,hy*2,nvar,nghost,naux,auxcrse)
          endif
      enddo

      do i  = nghost+1, mi2tot-nghost
        do j = nghost+1, mj2tot-nghost
          i_val = i-nghost
          j_val = j-nghost

!         Saving calculated errors to adjoint_module to be used 
!         for calculating tolerance in subsequent time steps
          errors(eptr(jg)+(i_val-1)*ny/2+j_val) = auxcrse(innerprod_index,i,j)

          if (auxcrse(innerprod_index,i,j) .ge. levtol(levm)) then
              rctcrse(1,i,j)  = DOFLAG
          else
              rctcrse(1,i,j)  = DONTFLAG
          endif

        enddo
      enddo

!     Updating fine grid inner product and flag values from coarse 
!     temporary vectors
      jfine   = nghost+1
      do j = nghost+1, mj2tot-nghost
        ifine   = nghost+1

        do i = nghost+1, mi2tot-nghost
          auxfine(innerprod_index,ifine,jfine) = auxcrse(innerprod_index,i,j)
          auxfine(innerprod_index,ifine+1,jfine) = auxcrse(innerprod_index,i,j)
          auxfine(innerprod_index,ifine,jfine+1) = auxcrse(innerprod_index,i,j)
          auxfine(innerprod_index,ifine+1,jfine+1) = auxcrse(innerprod_index,i,j)

          if (rctcrse(1,i,j) .eq. DOFLAG) then
!           ## never set rctflg to DONTFLAG, since flag2refine or
!           ## flagregions2 may have previously set it to DOFLAG
!           ## can only add DOFLAG pts in this routine
            rctflg(ifine,jfine)    = DOFLAG
            rctflg(ifine+1,jfine)  = DOFLAG
            rctflg(ifine,jfine+1)  = DOFLAG
            rctflg(ifine+1,jfine+1)= DOFLAG
          endif
          ifine   = ifine + 2
        enddo

        jfine   = jfine + 2
      enddo

    end subroutine errf1a

end module adjointsup_module
