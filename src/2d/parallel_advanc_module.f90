module parallel_advanc_module
    double precision :: dtcom,dxcom,dycom,tcom
    integer :: icom,jcom
contains
!  :::::::::::::: PAR_ADVANC :::::::::::::::::::::::::::::::::::::::::::
!  integrate this grid. grids are done in parallel.
!  extra subr. used to allow for stack based allocation of
!  flux arrays. They are only needed temporarily. If used alloc
!  array for them it has too long a lendim, makes too big
!  a checkpoint file, and is a big critical section.
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!> Integrate grid **mptr**. grids are done in parallel.
    subroutine par_advanc (mptr,mitot,mjtot,nvar,naux,dtnew)
        use amr_module
        use gauges_module, only: update_gauges, num_gauges
        implicit double precision (a-h,o-z)


        integer omp_get_thread_num, omp_get_max_threads
        integer mythread/0/, maxthreads/1/

        double precision fp(nvar,mitot,mjtot),fm(nvar,mitot,mjtot)
        double precision gp(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
        level = node(nestlevel,mptr)
        hx    = hxposs(level)
        hy    = hyposs(level)
        delt  = possk(level)
        nx    = node(ndihi,mptr) - node(ndilo,mptr) + 1
        ny    = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
        time  = rnode(timemult,mptr)

!$      mythread = omp_get_thread_num()

        locold = node(store2, mptr)
        locnew = node(store1, mptr)

!
!  copy old soln. values into  next time step's soln. values
!  since integrator will overwrite it. only for grids not at
!  the finest level. finest level grids do not maintain copies
!  of old and new time solution values.
!
        if (level .lt. mxnest) then
            ntot   = mitot * mjtot * nvar
!dir$ ivdep
            do 10 i = 1, ntot
 10             alloc(locold + i - 1) = alloc(locnew + i - 1)
        endif
!
        xlow = rnode(cornxlo,mptr) - nghost*hx
        ylow = rnode(cornylo,mptr) - nghost*hy

!$OMP CRITICAL(rv)
        rvol = rvol + nx * ny
        rvoll(level) = rvoll(level) + nx * ny
!$OMP END CRITICAL(rv)


        locaux = node(storeaux,mptr)
!
        if (node(ffluxptr,mptr) .ne. 0) then
            lenbc  = 2*(nx/intratx(level-1)+ny/intraty(level-1))
            locsvf = node(ffluxptr,mptr)
            locsvq = locsvf + nvar*lenbc
            locx1d = locsvq + nvar*lenbc
            call qad(alloc(locnew),mitot,mjtot,nvar, &
                     alloc(locsvf),alloc(locsvq),lenbc, &
                     intratx(level-1),intraty(level-1),hx,hy, &
                     naux,alloc(locaux),alloc(locx1d),delt,mptr)
        endif

!        # See if the grid about to be advanced has gauge data to output.
!        # This corresponds to previous time step, but output done
!        # now to make linear interpolation easier, since grid
!        # now has boundary conditions filled in.

!     should change the way print_gauges does io - right now is critical section
!     no more,  each gauge has own array.

        if (num_gauges > 0) then
            call update_gauges(alloc(locnew:locnew+nvar*mitot*mjtot), &
                               alloc(locaux:locaux+nvar*mitot*mjtot), &
                               xlow,ylow,nvar,mitot,mjtot,naux,mptr)
        endif

!
        if (dimensional_split .eq. 0) then
!           # Unsplit method
        ! call stepgrid(alloc(locnew),fm,fp,gm,gp, &
        !                 mitot,mjtot,nghost, &
        !                 delt,dtnew,hx,hy,nvar, &
        !                 xlow,ylow,time,mptr,naux,alloc(locaux))
        call stepgrid_cuda(alloc(locnew),fm,fp,gm,gp, &
                        mitot,mjtot,nghost, &
                        delt,dtnew,hx,hy,nvar, &
                        xlow,ylow,time,mptr,naux,alloc(locaux))
        else if (dimensional_split .eq. 1) then
!           # Godunov splitting
        call stepgrid_dimSplit(alloc(locnew),fm,fp,gm,gp, &
                     mitot,mjtot,nghost, &
                     delt,dtnew,hx,hy,nvar, &
                     xlow,ylow,time,mptr,naux,alloc(locaux))
        else 
!           # should never get here due to check in amr2
            write(6,*) '*** Strang splitting not supported'
            stop
        endif

        if (node(cfluxptr,mptr) .ne. 0) then
            call fluxsv(mptr,fm,fp,gm,gp, &
                     alloc(node(cfluxptr,mptr)),mitot,mjtot, &
                     nvar,listsp(level),delt,hx,hy)
        endif
        if (node(ffluxptr,mptr) .ne. 0) then
            lenbc = 2*(nx/intratx(level-1)+ny/intraty(level-1))
            locsvf = node(ffluxptr,mptr)
            call fluxad(fm,fp,gm,gp, &
                     alloc(locsvf),mptr,mitot,mjtot,nvar, &
                        lenbc,intratx(level-1),intraty(level-1), &
                     nghost,delt,hx,hy)
        endif
!
!        write(outunit,969) mythread,delt, dtnew
!969     format(" thread ",i4," updated by ",e15.7, " new dt ",e15.7)
        rnode(timemult,mptr)  = rnode(timemult,mptr)+delt
!
        return
    end subroutine par_advanc

! ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
! take a time step on a single grid. overwrite solution array q. 
! A modified version of the clawpack routine step2 is used.
!
! return fluxes in fm,fp and gm,gp.
! patch has room for ghost cells (mbc of them) around the grid.
! everything is the enlarged size (mitot by mjtot).
!
! mbc       = number of ghost cells  (= lwidth)
! mptr      = grid number  (for debugging)
! xlow,ylow = lower left corner of enlarged grid (including ghost cells).
! dt         = incoming time step
! dx,dy      = mesh widths for this grid
! dtnew      = return suggested new time step for this grid's soln.
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


    subroutine stepgrid(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dtnew,dx,dy, &
            nvar,xlow,ylow,time,mptr,maux,aux)

      use amr_module
      use cudafor
      ! use parallel_advanc_module, only: tcom
      implicit double precision (a-h,o-z)
      external rpn2,rpt2

      ! common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      parameter (msize=max1d+4)
      parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

      dimension q(nvar,mitot,mjtot)
      dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
      dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
! 	  dimension work(mwork)

      logical    debug,  dump
      data       debug/.false./,  dump/.false./


!
!     # set tcom = time.  This is in the common block comxyt that could
!     # be included in the Riemann solver, for example, if t is explicitly
!     # needed there.

      tcom = time

      if (dump) then
         write(outunit,*) "dumping grid ",mptr," at time ",time
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar) 
!    .                  ,(aux(ivar,i,j),ivar=1,maux)
 545        format(2i4,5e15.7)
         end do
         end do
      endif
!
      meqn   = nvar
      mx = mitot - 2*mbc
      my = mjtot - 2*mbc
      maxm = max(mx,my)       !# size for 1d scratch array
      mbig = maxm
      xlowmbc = xlow + mbc*dx
      ylowmbc = ylow + mbc*dy

!     # method(2:7) and mthlim
!     #    are set in the amr2ez file (read by amr)
!
      method(1) = 0

!
!
      call b4step2(mbc,mx,my,nvar,q, &
          xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)
!
!
!     # take one step on the conservation law:
!
      call step2(mbig,nvar,maux, &
          mbc,mx,my, &
          q,aux,dx,dy,dt,cflgrid, &
          fm,fp,gm,gp,rpn2,rpt2)
!
!
!       write(outunit,1001) mptr, node(nestlevel,mptr),cflgrid
!1001   format(' Courant # of grid', i4,
!    &        ' on level', i3, ' is  ', e10.3)
!

!$OMP  CRITICAL (cflm)

        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)

!
!       # update q
        dtdx = dt/dx
        dtdy = dt/dy
        do 50 j=mbc+1,mjtot-mbc
        do 50 i=mbc+1,mitot-mbc
        do 50 m=1,nvar
         if (mcapa.eq.0) then
!
!            # no capa array.  Standard flux differencing:

           q(m,i,j) = q(m,i,j) &
               - dtdx * (fm(m,i+1,j) - fp(m,i,j)) &
               - dtdy * (gm(m,i,j+1) - gp(m,i,j)) 
         else
!            # with capa array.
           q(m,i,j) = q(m,i,j) &
               - (dtdx * (fm(m,i+1,j) - fp(m,i,j)) &
               +  dtdy * (gm(m,i,j+1) - gp(m,i,j))) / aux(mcapa,i,j)
         endif

 50      continue
!
!
      if (method(5).eq.1) then
!        # with source term:   use Godunov splitting
         call src2(nvar,mbc,mx,my,xlowmbc,ylowmbc,dx,dy, &
             q,maux,aux,time,dt)
         endif
!
!
!
!     # output fluxes for debugging purposes:
      if (debug) then
         write(dbugunit,*)" fluxes for grid ",mptr
!        do 830 j = mbc+1, mjtot-1
            do 830 i = mbc+1, mitot-1
         do 830 j = mbc+1, mjtot-1
               write(dbugunit,831) i,j,fm(1,i,j),fp(1,i,j), &
                   gm(1,i,j),gp(1,i,j)
               do 830 m = 2, meqn
                  write(dbugunit,832) fm(m,i,j),fp(m,i,j), &
                      gm(m,i,j),gp(m,i,j)
  831          format(2i4,4d16.6)
  832          format(8x,4d16.6)
  830    continue
      endif

!
!
! For variable time stepping, use max speed seen on this grid to 
! choose the allowable new time step dtnew.  This will later be 
! compared to values seen on other grids.
!
       if (cflgrid .gt. 0.d0) then
           dtnew = dt*cfl/cflgrid
         else
!          # velocities are all zero on this grid so there's no 
!          # time step restriction coming from this grid.
            dtnew = rinfinity
          endif

!     # give a warning if Courant number too large...
!
      if (cflgrid .gt. cflv1) then
            write(*,810) cflgrid
            write(outunit,810) cflgrid, cflv1
  810       format('*** WARNING *** Courant number  =', d12.4, &
      '  is larger than input cfl_max = ', d12.4)
            endif
!
      if (dump) then
         write(outunit,*) "dumping grid ",mptr," after stepgrid"
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar)
         end do
         end do
      endif
      return
      end subroutine stepgrid

! ::::::::::::::::::: STEPGRID_CUDA ::::::::::::::::::::::::::::::::::::
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        subroutine stepgrid_cuda (q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dtnew,dx,dy, &
                    nvar,xlow,ylow,time,mptr,maux,aux)
        ! ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
        ! take a time step on a single grid. overwrite solution array q. 
        ! A modified version of the clawpack routine step2 is used.
        ! 
        ! return fluxes in fm,fp and gm,gp.
        ! patch has room for ghost cells (mbc of them) around the grid.
        ! everything is the enlarged size (mitot by mjtot).
        ! 
        ! mbc       = number of ghost cells  (= lwidth)
        ! mptr      = grid number  (for debugging)
        ! xlow,ylow = lower left corner of enlarged grid (including ghost cells).
        ! dt         = incoming time step
        ! dx,dy      = mesh widths for this grid
        ! dtnew      = return suggested new time step for this grid's soln.
        ! 
        ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

            use amr_module
            use cudafor
            implicit double precision (a-h,o-z)
            external rpn2,rpt2

            ! common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

            parameter (msize=max1d+4)
            parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

            dimension q(nvar,mitot,mjtot)
            dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
            dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
            dimension aux(maux,mitot,mjtot)

            ! variables for CUDA
            ! these also allocate space on device memory
            double precision, device :: q_d(nvar,mitot,mjtot)
            double precision, device :: fp_d(nvar,mitot,mjtot),gp_d(nvar,mitot,mjtot)
            double precision, device :: fm_d(nvar,mitot,mjtot),gm_d(nvar,mitot,mjtot)
            ! TODO: for now I assume maux = 0
            ! If I uncomment this statement below and maux = 0, the program has runtime 
            ! error since it cannot assign 0 bytes on device
            ! double precision, device :: aux_d(maux,mitot,mjtot)
            type(dim3) :: gridSize, blockSize 



            logical    debug,  dump
            data       debug/.false./,  dump/.false./

        ! 
        !     # set tcom = time.  This is in the common block comxyt that could
        !     # be included in the Riemann solver, for example, if t is explicitly
        !     # needed there.

            tcom = time

            ! if (dump) then
            !     write(outunit,*) "dumping grid ",mptr," at time ",time
            !     do i = 1, mitot
            !         do j = 1, mjtot
            !             write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar) 
            ! 545         format(2i4,5e15.7)
            !         end do
            !     end do
            ! endif
        ! 
            meqn   = nvar
            mx = mitot - 2*mbc
            my = mjtot - 2*mbc
            maxm = max(mx,my)       !# size for 1d scratch array
            mbig = maxm
            xlowmbc = xlow + mbc*dx
            ylowmbc = ylow + mbc*dy

            method(1) = 0
            ! TODO: write CUDA version of b4step2
            call b4step2(mbc,mx,my,nvar,q, &
                    xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)
            ! take one step on the conservation law:
            call step2(mbig,nvar,maux, &
                    mbc,mx,my, &
                    q,aux,dx,dy,dt,cflgrid, &
                    fm,fp,gm,gp,rpn2,rpt2)

            cfl_level = dmax1(cfl_level,cflgrid)


            dtdx = dt/dx
            dtdy = dt/dy

            ! copy data to device
            ! istat = cudaSetDevice(1)
            ! istat = cudaGetDevice( current_dev )
            ! write(*,*) 'We are currently using device: ',current_dev
            ! write(*,*) 'maux: ',maux
            istat = cudaMemcpy(q_d, q, nvar*mitot*mjtot, cudaMemcpyHostToDevice)
            istat = cudaMemcpy(fp_d, fp, nvar*mitot*mjtot, cudaMemcpyHostToDevice)
            istat = cudaMemcpy(fm_d, fm, nvar*mitot*mjtot, cudaMemcpyHostToDevice)
            istat = cudaMemcpy(gp_d, gp, nvar*mitot*mjtot, cudaMemcpyHostToDevice)
            istat = cudaMemcpy(gm_d, gm, nvar*mitot*mjtot, cudaMemcpyHostToDevice)
            ! istat = cudaMemcpy(aux_d, aux, maux*mitot*mjtot, cudaMemcpyHostToDevice)

            blockSize = dim3(16, 16, 1)
            gridSize = dim3(ceiling(real(mitot-2*mbc)/blockSize%x), &
                ceiling(real(mjtot-2*mbc)/blockSize%y), 1)
            call add_to_q<<<gridSize, blockSize>>> (q_d, fm_d, fp_d, gm_d, gp_d, &
                mbc, mcapa, mitot, mjtot, nvar, maux, &
                dtdx, dtdy)

            ! copy data back
            istat = cudaMemcpy(q, q_d, nvar*mitot*mjtot, cudaMemcpyDeviceToHost)
            istat = cudaMemcpy(fp, fp_d, nvar*mitot*mjtot, cudaMemcpyDeviceToHost)
            istat = cudaMemcpy(fm, fm_d, nvar*mitot*mjtot, cudaMemcpyDeviceToHost)
            istat = cudaMemcpy(gp, gp_d, nvar*mitot*mjtot, cudaMemcpyDeviceToHost)
            istat = cudaMemcpy(gm, gm_d, nvar*mitot*mjtot, cudaMemcpyDeviceToHost)
            ! istat = cudaMemcpy(aux, aux_d, maux*mitot*mjtot, cudaMemcpyDeviceToHost)


        ! 
        ! 
            if (method(5).eq.1) then
        !        # with source term:   use Godunov splitting
                call src2(nvar,mbc,mx,my,xlowmbc,ylowmbc,dx,dy, &
                        q,maux,aux,time,dt)
            endif
        ! 
        ! For variable time stepping, use max speed seen on this grid to 
        ! choose the allowable new time step dtnew.  This will later be 
        ! compared to values seen on other grids.
        ! 
            if (cflgrid .gt. 0.d0) then
                dtnew = dt*cfl/cflgrid
            else
        !          # velocities are all zero on this grid so there's no 
        !          # time step restriction coming from this grid.
                dtnew = rinfinity
            endif

            ! give a warning if Courant number too large...
            ! if (cflgrid .gt. cflv1) then
            !     write(*,810) cflgrid
            !     write(outunit,810) cflgrid, cflv1
            ! 810 format('*** WARNING *** Courant number  =', d12.4, &
            !     ' is larger than input cfl_max = ', d12.4)
            ! endif
        ! 
            return
        end subroutine stepgrid_cuda


        ! assume maux = 0
        attributes(global) subroutine add_to_q (q, fm, fp, gm, gp, &
                mbc, mcapa, mitot, mjtot, nvar, maux, &
                dtdx, dtdy)
            implicit none
            double precision, intent(out) :: q(nvar,mitot,mjtot)
            double precision, intent(in) :: fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
            double precision, intent(in) :: fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
            ! double precision, intent(in) :: aux(maux,mitot,mjtot)
            integer :: tx, ty, i, j, m
            ! TODO: maybe able to compute size of these arrays in this subroutine
            ! and thus avoid passing so many parameters
            ! passed by value
            integer, value :: mbc, mcapa, mitot, mjtot, nvar, maux
            double precision, value :: dtdx, dtdy
            tx = threadIdx%x 
            ty = threadIdx%y 
            i =  (blockIdx%x-1)*blockDim%x + tx + mbc
            j =  (blockIdx%y-1)*blockDim%y + ty + mbc
            ! each thread will handle all m unknowns of of a cell
            if ( i <= mitot-mbc .and. j <= mjtot-mbc ) then
                if (mcapa.eq.0) then
                ! no capa array.  Standard flux differencing:
                    do m = 1, nvar
                        q(m,i,j) = q(m,i,j) &
                            - dtdx * (fm(m,i+1,j) - fp(m,i,j)) &
                            - dtdy * (gm(m,i,j+1) - gp(m,i,j)) 
                    end do
                else
                ! with capa array.
                    write(*,*) "Not handling mcapa > 0"
                    ! do m = 1, nvar
                    !     q(m,i,j) = q(m,i,j) &
                    !         - (dtdx * (fm(m,i+1,j) - fp(m,i,j)) &
                    !         +  dtdy * (gm(m,i,j+1) - gp(m,i,j))) / aux(mcapa,i,j)
                    ! end do
                endif
            endif
        end subroutine add_to_q
end module parallel_advanc_module
