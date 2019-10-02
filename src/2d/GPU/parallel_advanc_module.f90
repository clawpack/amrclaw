module parallel_advanc_module

    use amr_module

    real(CLAW_REAL) :: dtcom,dxcom,dycom,tcom
    integer :: icom,jcom
contains

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

!
! This subroutine should be used only if CUDA is defined
!
    subroutine stepgrid_soa(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dx,dy, &
            nvar,xlow,ylow,time,mptr,maux,aux,ngrids,id,cfls)

        use amr_module
#ifdef CUDA
        use memory_module, only: gpu_allocate, gpu_deallocate, cpu_allocate_pinned, cpu_deallocate_pinned
        use cuda_module, only: device_id, wait_for_all_gpu_tasks
        use cuda_module, only: get_cuda_stream
        use step2_cuda_module, only: step2_fused
        use cudafor
#endif
        implicit real(CLAW_REAL) (a-h,o-z)

        parameter (msize=max1d+4)
        parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

        integer, intent(in) :: id, ngrids
        real(CLAW_REAL) ::  xlow, ylow
        real(CLAW_REAL), intent(in) :: dt,dx,dy
        real(CLAW_REAL), intent(out) :: cfls(SPACEDIM,ngrids)
        ! These are all in SoA format
        real(CLAW_REAL) ::   q(mitot,mjtot,nvar)
        real(CLAW_REAL) ::  fp(mitot,mjtot,nvar),gp(mitot,mjtot,nvar)
        real(CLAW_REAL) ::  fm(mitot,mjtot,nvar),gm(mitot,mjtot,nvar)
        real(CLAW_REAL) :: aux(mitot,mjtot,maux)
#ifdef CUDA
        attributes(device) :: q
        attributes(device) :: fp, fm, gp, gm
        attributes(device) :: cfls
        ! attributes(device) :: aux
#endif

        real(CLAW_REAL) :: dtdx, dtdy
        integer :: i,j,m

        logical    debug,  dump
        data       debug/.false./,  dump/.false./

!     # set tcom = time.  This is in the common block comxyt that could
!     # be included in the Riemann solver, for example, if t is explicitly
!     # needed there.

        tcom = time

#ifndef CUDA
        if (dump) then
           write(outunit,*) "dumping grid ",mptr," at time ",time
           do i = 1, mitot
           do j = 1, mjtot
              write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar) 
!    .                    ,(aux(ivar,i,j),ivar=1,maux)
 545          format(2i4,5e15.7)
           end do
           end do
        endif
#endif
!
        meqn   = nvar
        mx = mitot - 2*mbc
        my = mjtot - 2*mbc
        maxm = max(mx,my)       !# size for 1d scratch array
        mbig = maxm
        xlowmbc = xlow + mbc*dx
        ylowmbc = ylow + mbc*dy

!       # method(2:7) and mthlim
!       #    are set in the amr2ez file (read by amr)
!
        method(1) = 0

!
!
#ifndef CUDA
        ! For now this does nothing
        call b4step2(mbc,mx,my,nvar,q, &
          xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)
#endif

! assume no aux here
      call step2_fused(mbig,nvar,maux, &
          mbc,mx,my, &
          q,dx,dy,dt,cfls, &
          fm,fp,gm,gp,mptr,ngrids,id)

!       # update q
        dtdx = dt/dx
        dtdy = dt/dy

        !$cuf kernel do(3) <<<*, *, 0, get_cuda_stream(id,device_id)>>>
        do m=1,nvar
            do j=mbc+1,mjtot-mbc
                do i=mbc+1,mitot-mbc
                    if (mcapa.eq.0) then
                        !            # no capa array.  Standard flux differencing:
                        q(i,j,m) = q(i,j,m) &
                            - dtdx * (fm(i+1,j,m) - fp(i,j,m)) &
                            - dtdy * (gm(i,j+1,m) - gp(i,j,m)) 
                    else
                        print *, "With-capa-array case is not implemented"
                        stop
                        !            # with capa array.
                        ! q(m,i,j) = q(m,i,j) &
                        !     - (dtdx * (fm(m,i+1,j) - fp(m,i,j)) &
                        !     +  dtdy * (gm(m,i,j+1) - gp(m,i,j))) / aux(mcapa,i,j)
                    endif
                enddo
            enddo
        enddo
!
!
#ifndef CUDA
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
!           do 830 j = mbc+1, mjtot-1
               do 830 i = mbc+1, mitot-1
            do 830 j = mbc+1, mjtot-1
                  write(dbugunit,831) i,j,fm(1,i,j),fp(1,i,j), &
                      gm(1,i,j),gp(1,i,j)
                  do 830 m = 2, meqn
                     write(dbugunit,832) fm(m,i,j),fp(m,i,j), &
                         gm(m,i,j),gp(m,i,j)
  831             format(2i4,4d16.6)
  832             format(8x,4d16.6)
  830       continue
        endif
#endif

        return
    end subroutine stepgrid_soa

!
! This subroutine should be used only if CUDA is defined
!
#ifdef CUDA
    subroutine stepgrid_dimsplit_soa(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dx,dy, &
            nvar,xlow,ylow,time,mptr,maux,aux,ngrids,id,cfls)

        use amr_module
        use memory_module, only: gpu_allocate, gpu_deallocate, cpu_allocate_pinned, cpu_deallocate_pinned
        use cuda_module, only: device_id, wait_for_all_gpu_tasks
        use cuda_module, only: get_cuda_stream
        use step2_cuda_module, only: step2_and_update
        use cudafor
        implicit real(CLAW_REAL) (a-h,o-z)

        parameter (msize=max1d+4)
        parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

        integer, intent(in) :: id, ngrids
        real(CLAW_REAL) ::  xlow, ylow
        real(CLAW_REAL), intent(in) :: dt,dx,dy
        real(CLAW_REAL), intent(out) :: cfls(SPACEDIM,ngrids)
        ! These are all in SoA format
        real(CLAW_REAL) ::   q(mitot,mjtot,nvar)
        real(CLAW_REAL) ::  fp(mitot,mjtot,nvar),gp(mitot,mjtot,nvar)
        real(CLAW_REAL) ::  fm(mitot,mjtot,nvar),gm(mitot,mjtot,nvar)
        real(CLAW_REAL) :: aux(mitot,mjtot,maux)
        attributes(device) :: q
        attributes(device) :: fp, fm, gp, gm
        attributes(device) :: cfls
        ! attributes(device) :: aux

        real(CLAW_REAL) :: dtdx, dtdy
        integer :: i,j,m

        logical    debug,  dump
        data       debug/.false./,  dump/.false./

!     # set tcom = time.  This is in the common block comxyt that could
!     # be included in the Riemann solver, for example, if t is explicitly
!     # needed there.

        meqn   = nvar
        mx = mitot - 2*mbc
        my = mjtot - 2*mbc
        maxm = max(mx,my)       !# size for 1d scratch array
        mbig = maxm
        xlowmbc = xlow + mbc*dx
        ylowmbc = ylow + mbc*dy

!       # method(2:7) and mthlim
!       #    are set in the amr2ez file (read by amr)
!
        method(1) = 0

!
!
        ! For now this does nothing
        call b4step2(mbc,mx,my,nvar,q, &
          xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)

! assume no aux here
      call step2_and_update(mbig,nvar,maux, &
          mbc,mx,my, &
          q,dx,dy,dt,cfls, &
          fm,fp,gm,gp,mptr,ngrids,id)

!       # update q
        dtdx = dt/dx
        dtdy = dt/dy

        !$cuf kernel do(3) <<<*, *, 0, get_cuda_stream(id,device_id)>>>
        do m=1,nvar
            do j=mbc+1,mjtot-mbc
                do i=mbc+1,mitot-mbc
                    if (mcapa.eq.0) then
                        !            # no capa array.  Standard flux differencing:
                        q(i,j,m) = q(i,j,m) &
                            - dtdx * (fm(i+1,j,m) - fp(i,j,m)) &
                            - dtdy * (gm(i,j+1,m) - gp(i,j,m)) 
                    else
                        print *, "With-capa-array case is not implemented"
                        stop
                        !            # with capa array.
                        ! q(m,i,j) = q(m,i,j) &
                        !     - (dtdx * (fm(m,i+1,j) - fp(m,i,j)) &
                        !     +  dtdy * (gm(m,i,j+1) - gp(m,i,j))) / aux(mcapa,i,j)
                    endif
                enddo
            enddo
        enddo
!
!
        ! this has not been implemented for CUDA
        ! if (method(5).eq.1) then
!       !  # with source term:   use Godunov splitting
        !  call src2(nvar,mbc,mx,my,xlowmbc,ylowmbc,dx,dy, &
        !      q,maux,aux,time,dt)
        ! endif
!
!
!
!         ! this has not been implemented for CUDA
! !     # output fluxes for debugging purposes:
!         if (debug) then
!             write(dbugunit,*)" fluxes for grid ",mptr
! !           do 830 j = mbc+1, mjtot-1
!                do 830 i = mbc+1, mitot-1
!             do 830 j = mbc+1, mjtot-1
!                   write(dbugunit,831) i,j,fm(1,i,j),fp(1,i,j), &
!                       gm(1,i,j),gp(1,i,j)
!                   do 830 m = 2, meqn
!                      write(dbugunit,832) fm(m,i,j),fp(m,i,j), &
!                          gm(m,i,j),gp(m,i,j)
!   831             format(2i4,4d16.6)
!   832             format(8x,4d16.6)
!   830       continue
!         endif

        return
    end subroutine stepgrid_dimsplit_soa

#endif
end module parallel_advanc_module
