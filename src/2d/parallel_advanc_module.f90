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
        call stepgrid(alloc(locnew),fm,fp,gm,gp, &
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
#ifdef CUDA
        use memory_module, only: gpu_allocate, gpu_deallocate, cpu_allocate_pinned, cpu_deallocated_pinned
        use cuda_module, only: device_id, wait_for_all_gpu_tasks, cuda_streams
        use cudafor
#endif
        implicit double precision (a-h,o-z)
        external rpn2,rpt2

        parameter (msize=max1d+4)
        parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

        dimension q(nvar,mitot,mjtot)
        dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
        dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
        dimension aux(maux,mitot,mjtot)

        double precision :: dtdx, dtdy
        integer :: i,j,m

#ifdef CUDA
        ! These are all in SoA format
        double precision, dimension(:,:,:), pointer, contiguous :: &
            q1, fp1, gp1, fm1, gm1
        double precision, dimension(:,:,:), pointer, contiguous, device :: &
            q_d, fp_d, gp_d, fm_d, gm_d
        integer :: data_size, aux_size
        integer :: cudaResult
#endif
        logical    debug,  dump
        data       debug/.false./,  dump/.false./

#ifdef CUDA
        call cpu_allocate_pinned( q1, 1, mitot, 1, mjtot, 1, nvar) 
        call cpu_allocate_pinned(fp1, 1, mitot, 1, mjtot, 1, nvar) 
        call cpu_allocate_pinned(gp1, 1, mitot, 1, mjtot, 1, nvar) 
        call cpu_allocate_pinned(fm1, 1, mitot, 1, mjtot, 1, nvar) 
        call cpu_allocate_pinned(gm1, 1, mitot, 1, mjtot, 1, nvar) 

        call gpu_allocate( q_d, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(fp_d, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(gp_d, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(fm_d, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        call gpu_allocate(gm_d, device_id, 1, mitot, 1, mjtot, 1, nvar) 
        ! if (maux > 0) then
        !     call gpu_allocate(aux_d, device_id, 1, maux, 1, mitot, 1, mjtot) 
        ! endif
        data_size = nvar*mitot*mjtot 
        ! aux_size = maux*mitot*mjtot 
#endif

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

#ifdef CUDA
        ! convert q to SoA
        do i = 1, mitot
            do j = 1,mjtot
                do m = 1,nvar
                    q1(i,j,m) = q(m,i,j)
                enddo
            enddo
        enddo
        cudaResult = cudaMemcpyAsync(q_d,  q1, data_size, cudaMemcpyHostToDevice, cuda_streams(1,device_id))
        ! if (maux > 0) then
        !     cudaResult = cudaMemcpy(aux_d,aux, aux_size, cudaMemcpyHostToDevice)
        ! endif
#endif

! assume no aux here
#ifdef CUDA
      call step2_fused(mbig,nvar,maux, &
          mbc,mx,my, &
          q_d,dx,dy,dt,cflgrid, &
          fm_d,fp_d,gm_d,gp_d,rpn2,rpt2)
#else
      call step2_fused(mbig,nvar,maux, &
          mbc,mx,my, &
          q,dx,dy,dt,cflgrid, &
          fm,fp,gm,gp,rpn2,rpt2)
#endif

!$OMP  CRITICAL (cflm)

        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)

!
!       # update q
        dtdx = dt/dx
        dtdy = dt/dy

#ifdef CUDA
        !$cuf kernel do(3) <<<*, *>>>
        do m=1,nvar
            do j=mbc+1,mjtot-mbc
                do i=mbc+1,mitot-mbc
                    if (mcapa.eq.0) then
                        !            # no capa array.  Standard flux differencing:
                        q_d(i,j,m) = q_d(i,j,m) &
                            - dtdx * (fm_d(i+1,j,m) - fp_d(i,j,m)) &
                            - dtdy * (gm_d(i,j+1,m) - gp_d(i,j,m)) 
                    else
                        print *, "With-capa-array case is not implemented"
                        !            # with capa array.
                        ! q_d(m,i,j) = q_d(m,i,j) &
                        !     - (dtdx * (fm_d(m,i+1,j) - fp_d(m,i,j)) &
                        !     +  dtdy * (gm_d(m,i,j+1) - gp_d(m,i,j))) / aux_d(mcapa,i,j)
                    endif
                enddo
            enddo
        enddo
#else
        do j=mbc+1,mjtot-mbc
            do i=mbc+1,mitot-mbc
                do m=1,nvar
                    if (mcapa.eq.0) then
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
                enddo
            enddo
        enddo
#endif

#ifdef CUDA
        cudaResult = cudaMemcpyAsync( q1, q_d,data_size, cudaMemcpyDeviceToHost, cuda_streams(1,device_id))
        cudaResult = cudaMemcpyAsync(fm1,fm_d,data_size, cudaMemcpyDeviceToHost, cuda_streams(1,device_id))
        cudaResult = cudaMemcpyAsync(fp1,fp_d,data_size, cudaMemcpyDeviceToHost, cuda_streams(1,device_id))
        cudaResult = cudaMemcpyAsync(gm1,gm_d,data_size, cudaMemcpyDeviceToHost, cuda_streams(1,device_id))
        cudaResult = cudaMemcpyAsync(gp1,gp_d,data_size, cudaMemcpyDeviceToHost, cuda_streams(1,device_id))
        ! if (maux > 0) then
        !     cudaResult = cudaMemcpy(aux,aux_d,aux_size, cudaMemcpyDeviceToHost)
        ! endif

        call wait_for_all_gpu_tasks(device_id)
#endif

#ifdef CUDA
        ! convert q1 back to AoS
        do i = 1, mitot
            do j = 1,mjtot
                do m = 1,nvar
                    q(m,i,j) = q1(i,j,m)
                enddo
            enddo
        enddo
        ! convert fm, fp, gm, gp to AoS
        do i = 1, mitot
            do j = 1,mjtot
                do m = 1,nvar
                    fm(m,i,j) = fm1(i,j,m)
                    fp(m,i,j) = fp1(i,j,m)
                    gm(m,i,j) = gm1(i,j,m)
                    gp(m,i,j) = gp1(i,j,m)
                enddo
            enddo
        enddo

        call cpu_deallocated_pinned( q1) 
        call cpu_deallocated_pinned(fp1) 
        call cpu_deallocated_pinned(gp1) 
        call cpu_deallocated_pinned(fm1) 
        call cpu_deallocated_pinned(gm1) 

        call gpu_deallocate( q_d, device_id) 
        call gpu_deallocate(fp_d, device_id) 
        call gpu_deallocate(gp_d, device_id) 
        call gpu_deallocate(fm_d, device_id) 
        call gpu_deallocate(gm_d, device_id) 
        ! if (maux > 0) then
        !     call gpu_deallocate(aux_d) 
        ! endif
#endif

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

end module parallel_advanc_module
