module stepgrid_module
    use amr_module
#ifdef CUDA
    interface
        subroutine stepgrid_cudaclaw(cellsX, cellsY, ghostCells, &
            startX, endX, startY, endY, dt, &
            q_tmp, qNew, &
            coefficients, &
            waveSpeedsX, waveSpeedsY, &
            numStates, numCoefficients, &
            cfls_d, ngrids, mcapa, id, dev_id) &
            bind(C,NAME='call_C_limited_riemann_update')


            use iso_c_binding

            implicit none
            integer(kind=c_int), value, intent(in) :: cellsX, cellsY, ghostCells, &
                ngrids, id, dev_id, numStates, numCoefficients, mcapa

            #if (CLAW_REAL == 8)
                real(kind=c_double), value, intent(in) :: startX, endX, startY, endY, dt
                real(kind=c_double), device :: q_tmp(cellsX,cellsY,numStates)
                real(kind=c_double), device :: qNew(cellsX,cellsY,numStates)
                real(kind=c_double), device :: coefficients(cellsX,cellsY,numStates)
                real(kind=c_double), device :: cfls_d(SPACEDIM)
                real(kind=c_double), device :: waveSpeedsX(ws_len)
                real(kind=c_double), device :: waveSpeedsY(ws_len)
            #else
                real(kind=c_float), value, intent(in) :: startX, endX, startY, endY, dt
                real(kind=c_float), device :: q_tmp(cellsX,cellsY,numStates)
                real(kind=c_float), device :: qNew(cellsX,cellsY,numStates)
                real(kind=c_float), device :: coefficients(cellsX,cellsY,numStates)
                real(kind=c_float), device :: cfls_d(SPACEDIM)
                real(kind=c_float), device :: waveSpeedsX(ws_len)
                real(kind=c_float), device :: waveSpeedsY(ws_len)
            #endif

        end subroutine stepgrid_cudaclaw
    end interface
#endif
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
#ifndef CUDA
    subroutine par_advanc (mptr,mitot,mjtot,nvar,naux,dtnew)
        use amr_module
        use gauges_module, only: update_gauges, num_gauges
        implicit real(CLAW_REAL) (a-h,o-z)


        integer omp_get_thread_num, omp_get_max_threads
        integer mythread/0/, maxthreads/1/

        real(CLAW_REAL) fp(nvar,mitot,mjtot),fm(nvar,mitot,mjtot)
        real(CLAW_REAL) gp(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
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
#endif
end module stepgrid_module
