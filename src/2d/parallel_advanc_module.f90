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
        use memory_module, only: gpu_allocate, gpu_deallocate
        use cuda_module, only: device_id, wait_for_all_gpu_tasks
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
        double precision, dimension(:,:,:), pointer, contiguous, device :: &
            q_d, fp_d, gp_d, fm_d, gm_d, aux_d
        integer :: data_size, aux_size
        integer :: cudaResult
#endif
        ! 	  dimension work(mwork)

        logical    debug,  dump
        data       debug/.false./,  dump/.false./

#ifdef CUDA
        call gpu_allocate(q_d, device_id, 1, nvar, 1, mitot, 1, mjtot) 
        call gpu_allocate(fp_d, device_id, 1, nvar, 1, mitot, 1, mjtot) 
        call gpu_allocate(gp_d, device_id, 1, nvar, 1, mitot, 1, mjtot) 
        call gpu_allocate(fm_d, device_id, 1, nvar, 1, mitot, 1, mjtot) 
        call gpu_allocate(gm_d, device_id, 1, nvar, 1, mitot, 1, mjtot) 
        call gpu_allocate(aux_d, device_id, 1, maux, 1, mitot, 1, mjtot) 
        data_size = nvar*mitot*mjtot 
        aux_size = maux*mitot*mjtot 
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

#ifdef CUDA
        cudaResult = cudaMemcpy(  q_d,  q,data_size, cudaMemcpyHostToDevice)
        cudaResult = cudaMemcpy( fm_d, fm,data_size, cudaMemcpyHostToDevice)
        cudaResult = cudaMemcpy( fp_d, fp,data_size, cudaMemcpyHostToDevice)
        cudaResult = cudaMemcpy( gm_d, gm,data_size, cudaMemcpyHostToDevice)
        cudaResult = cudaMemcpy( gp_d, gp,data_size, cudaMemcpyHostToDevice)
        cudaResult = cudaMemcpy(aux_d,aux, aux_size, cudaMemcpyHostToDevice)
#endif
!
!       # update q
        dtdx = dt/dx
        dtdy = dt/dy

#ifdef CUDA
        !$cuf kernel do(3) <<<*, *>>>
        do j=mbc+1,mjtot-mbc
            do i=mbc+1,mitot-mbc
                do m=1,nvar
                    if (mcapa.eq.0) then
                        !            # no capa array.  Standard flux differencing:
                        q_d(m,i,j) = q_d(m,i,j) &
                            - dtdx * (fm_d(m,i+1,j) - fp_d(m,i,j)) &
                            - dtdy * (gm_d(m,i,j+1) - gp_d(m,i,j)) 
                    else
                        !            # with capa array.
                        q_d(m,i,j) = q_d(m,i,j) &
                            - (dtdx * (fm_d(m,i+1,j) - fp_d(m,i,j)) &
                            +  dtdy * (gm_d(m,i,j+1) - gp_d(m,i,j))) / aux_d(mcapa,i,j)
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
        cudaResult = cudaMemcpy(q,q_d,  data_size, cudaMemcpyDeviceToHost)
        cudaResult = cudaMemcpy(fm,fm_d,data_size, cudaMemcpyDeviceToHost)
        cudaResult = cudaMemcpy(fp,fp_d,data_size, cudaMemcpyDeviceToHost)
        cudaResult = cudaMemcpy(gm,gm_d,data_size, cudaMemcpyDeviceToHost)
        cudaResult = cudaMemcpy(gp,gp_d,data_size, cudaMemcpyDeviceToHost)
        cudaResult = cudaMemcpy(aux,aux_d,aux_size, cudaMemcpyDeviceToHost)
#endif

#ifdef CUDA
        call gpu_deallocate(q_d) 
        call gpu_deallocate(fp_d) 
        call gpu_deallocate(gp_d) 
        call gpu_deallocate(fm_d) 
        call gpu_deallocate(gm_d) 
        call gpu_deallocate(aux_d) 
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
