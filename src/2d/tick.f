c
c  -------------------------------------------------------------
c
      subroutine tick(nvar,cut,nstart,vtime,time,naux,start_time,
     &                rest,dt_max)
c
      use amr_module
      use gauges_module, only: setbestsrc, num_gauges
      use gauges_module, only: print_gauges_and_reset_nextLoc 

      implicit double precision (a-h,o-z)
c     include  "call.i"

      logical    vtime, dumpout/.false./, dumpchk/.false./
      logical    rest, dump_final, stopFound
      dimension dtnew(maxlv), ntogo(maxlv), tlevel(maxlv)
      integer(kind=8) ::   clock_start, clock_finish, clock_rate
      integer(kind=8) ::   tick_clock_finish, tick_clock_rate  
      real(kind=8) :: cpu_start,cpu_finish

c
c :::::::::::::::::::::::::::: TICK :::::::::::::::::::::::::::::
!> main driver routine.  controls:
!!       integration  of all grids.
!!       error estimation / regridding
!!       output counting
!!       updating of fine to coarse grids

c  parameters:
c     nstop   = # of coarse grid time steps to be taken
c     iout    = output interval every 'iout' coarse time steps
c               (if 0, not used - set to inf.)
c     vtime   = true for variable timestep, calculated each coarse step
c
c  integration strategy is to advance a fine grid until it catches
c  up to the coarse grid. this strategy is applied recursively.
c  coarse grid goes first.
c
c  nsteps: used to count how number steps left for a level to be
c          integrated before it catches up with the next coarser level.
c  ncycle: counts number of coarse grid steps = # cycles.
c
c  icheck: counts the number of steps (incrementing by 1
c          each step) to keep track of when that level should
c          have its error estimated and finer levels should be regridded.
c ::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::
c


      ncycle         = nstart
      call setbestsrc()     ! need at very start of run, including restart
      if (iout .eq. 0) then
c        # output_style 1 or 2
         iout  = iinfinity
         nextout = 0
         if (nout .gt. 0) then
            nextout = 1
            if (nstart .gt. 0) then
c              # restart: make sure output times start after restart time
               do ii = 1, nout
                 if (tout(ii) .gt. time) then
                   nextout = ii
                   go to 2
                 endif
               end do
  2         continue
            endif
         endif
      endif

      nextchk = 1
      if ((nstart .gt. 0) .and. (abs(checkpt_style).eq.2)) then
c        if this is a restart, make sure chkpt times start after restart time
         do ii = 1, nchkpt
           if (tchk(ii) .gt. time) then
              nextchk = ii
              go to 3
              endif
           enddo
  3      continue
         endif

      tlevel(1)      = time

      do 5 i       = 2, mxnest
       tlevel(i) = tlevel(1)
 5     continue
c
c  ------ start of coarse grid integration loop. ------------------
c
 20   if (ncycle .ge. nstop .or. time .ge. tfinal) goto 999

      if (nout .gt. 0) then
          if (nextout  .le. nout) then
             outtime       = tout(nextout)
          else
             outtime       = rinfinity
          endif
      else
          outtime = tfinal
      endif

      if (nextchk  .le. nchkpt) then
         chktime       = tchk(nextchk)
      else
         chktime       = rinfinity
      endif

      dumpout = .false.  !# may be reset below

      if (time.lt.outtime .and. time+1.001*possk(1) .ge. outtime) then
c        ## adjust time step  to hit outtime exactly, and make output
c        #  apr 2010 mjb: modified so allow slightly larger timestep to
c        #  hit output time exactly, instead of taking minuscule timestep
c        #  should still be stable since increase dt in only 3rd digit.
         oldposs = possk(1)
         possk(1) = outtime - time
c        write(*,*)" old possk is ", possk(1)
         diffdt = oldposs - possk(1)  ! if positive new step is smaller


         if (.false.) then  
            write(*,122) diffdt,outtime  ! notify of change
 122        format(" Adjusting timestep by ",e10.3,
     .             " to hit output time of ",e13.6)
c           write(*,*)" new possk is ", possk(1)
            if (diffdt .lt. 0.) then ! new step is slightly larger
              pctIncrease = -100.*diffdt/oldposs   ! minus sign to make whole expr. positive
              write(*,123) pctIncrease
 123          format(" New step is ",f8.2," % larger.",
     .               "  Should still be stable")
              endif
            endif


         do i = 2, mxnest
            possk(i) = possk(i-1) / kratio(i-1)
            enddo
         if (nout .gt. 0) then
            nextout = nextout + 1
            dumpout = .true.
            endif
      endif


      if (time.lt.chktime .and. time + possk(1) .ge. chktime) then
c        ## adjust time step  to hit chktime exactly, and do checkpointing
         possk(1) = chktime - time
         do 13 i = 2, mxnest
            possk(i) = possk(i-1) / kratio(i-1)
 13      continue
         nextchk = nextchk + 1
        dumpchk = .true.
      else
        dumpchk = .false.
      endif

c
      level        = 1
      ntogo(level) = 1
      do i = 1, maxlv
         dtnew(i)  = rinfinity
      enddo
c
c     ------------- regridding  time?  ---------
c
c check if either
c   (i)  this level should have its error estimated before being advanced
c   (ii) this level needs to provide boundary values for either of
c        next 2 finer levels to have their error estimated.
c        this only affects two grid levels higher, occurs because
c        previous time step needs boundary vals for giant step.
c  no error estimation on finest possible grid level
c
 60       continue
          if (icheck(level) .ge. kcheck) then
               lbase = level
          else if (level+1 .ge. mxnest) then
               go to 90
          else if (icheck(level+1) .ge. kcheck) then
               lbase = level+1
          else if (level+2 .ge. mxnest) then
               go to 90
          else if (icheck(level+2) .ge. kcheck) then
               lbase = level+2
          else
               go to 90
          endif
          if (lbase .eq. mxnest .or. lbase .gt. lfine) go to 70
c
c regrid level 'lbase+1' up to finest level.
c level 'lbase' stays fixed.
c
          if (rprint) write(outunit,101) lbase
101       format(8h  level ,i5,32h  stays fixed during regridding )

          call system_clock(clock_start,clock_rate)
          call cpu_time(cpu_start)
          call regrid(nvar,lbase,cut,naux,start_time)
          call system_clock(clock_finish,clock_rate)
          call cpu_time(cpu_finish)
          timeRegridding = timeRegridding + clock_finish - clock_start
          timeRegriddingCPU=timeRegriddingCPU+cpu_finish-cpu_start

          call setbestsrc()     ! need at every grid change
c         call outtre(lstart(lbase+1),.false.,nvar,naux)
c note negative time to signal regridding output in plots
c         call valout(lbase,lfine,-tlevel(lbase),nvar,naux)
c
c  maybe finest level in existence has changed. reset counters.
c
          if (rprint .and. lbase .lt. lfine) then
             call outtre(lstart(lbase+1),.false.,nvar,naux)
          endif
 70       continue
          do 80  i  = lbase, lfine
             icheck(i) = 0
 80       continue
          do 81  i  = lbase+1, lfine
             tlevel(i) = tlevel(lbase)
 81       continue
c
c          MJB: modified to check level where new grids start, which is lbase+1
          !if (verbosity_regrid.ge.lbase+1) then
          if (.false.) then  ! don't need to print these every time:
                 do levnew = lbase+1,lfine
                     write(6,1006) intratx(levnew-1),intraty(levnew-1),
     &                             kratio(levnew-1),levnew
 1006                format('   Refinement ratios...  in x:', i3, 
     &                 '  in y:',i3,'  in t:',i3,' for level ',i4)
                 end do

              endif

c  ------- done regridding --------------------
c
c integrate all grids at level 'level'.
c
 90       continue


          call advanc(level,nvar,dtlevnew,vtime,naux)

c         # rjl modified 6/17/05 to print out *after* advanc and print cfl
c         # rjl & mjb changed to cfl_level, 3/17/10

          timenew = tlevel(level)+possk(level)
          if (tprint) then
              write(outunit,100)level,cfl_level,possk(level),timenew
              endif
          if (method(4).ge.level) then
              write(6,100)level,cfl_level,possk(level),timenew
              endif
100       format(' AMRCLAW: level ',i2,'  CFL = ',e10.3,
     &           '  dt = ',e11.4,  '  final t = ',e13.6)


c        # to debug individual grid updates...
c        call valout(level,level,time,nvar,naux)
c
c done with a level of integration. update counts, decide who next.
c
          ntogo(level)  = ntogo(level) - 1
          dtnew(level)  = dmin1(dtnew(level),dtlevnew)
          tlevel(level) = tlevel(level) + possk(level)
          icheck(level) = icheck(level) + 1
c
          if (level .lt. lfine) then
             level = level + 1
c            #  check if should adjust finer grid time step to start wtih
             if (((possk(level-1) - dtnew(level-1))/dtnew(level-1)) .gt.
     .            .05) then
                dttemp = dtnew(level-1)/kratio(level-1)
                ntogo(level) = (tlevel(level-1)-tlevel(level))/dttemp+.9
              else
                ntogo(level) = kratio(level-1)
              endif
             possk(level) = possk(level-1)/ntogo(level)
             go to 60
          endif
c
 105      if (level .eq. 1) go to 110
              if (ntogo(level) .gt. 0) then
c                same level goes again. check for ok time step
 106             if ((possk(level)-dtnew(level))/dtnew(level)
     .                .gt. .05)  then
c                   adjust time steps for this and finer levels
                    ntogo(level) = ntogo(level) + 1
                    possk(level) = (tlevel(level-1)-tlevel(level))/
     .                             ntogo(level)
                    go to 106
                 endif
                 go to 60
              else
                 level = level - 1
                 call system_clock(clock_start,clock_rate)
                 call update(level,nvar,naux)
                 call system_clock(clock_finish,clock_rate)
                 timeUpdating=timeUpdating+clock_finish-clock_start
              endif
          go to 105
c
c  --------------one complete coarse grid integration cycle done. -----
c
c      time for output?  done with the whole thing?
c
 110      continue
          time    = time   + possk(1)
          ncycle  = ncycle + 1
          call conck(1,nvar,naux,time,rest)

      if ( vtime) then
c
c         find new dt for next cycle (passed back from integration routine).
           do 115 i = 2, lfine
             ii = lfine+1-i
             dtnew(ii) = min(dtnew(ii),dtnew(ii+1)*kratio(ii))
 115       continue
c          make sure not to exceed largest permissible dt
           dtnew(1) = min(dtnew(1),dt_max)  
           possk(1) = dtnew(1)    ! propagate new timestep to hierarchy
           do 120 i = 2, mxnest
             possk(i) = possk(i-1) / kratio(i-1)
 120       continue

      endif

       if ((abs(checkpt_style).eq.3 .and. 
     &      mod(ncycle,checkpt_interval).eq.0) .or. dumpchk) then
                call check(ncycle,time,nvar,naux)
                dumpchk = .true.
       endif

       if ((mod(ncycle,iout).eq.0) .or. dumpout) then
         call valout(1,lfine,time,nvar,naux)
         if (abs(checkpt_style).eq.4) then
            call check(ncycle,time,nvar,naux)
            dumpchk = .true.
         endif
         if (printout) call outtre(mstart,.true.,nvar,naux)
       endif

       ! new STOP feature to do immediate checkpt and exit
       inquire(FILE="STOP.txt",exist=stopFound) 
          if (stopFound) then      
          write(*,*)"STOP file found. Checkpointing and Stopping"
          write(*,*)"REMEMBER to remove file before restarting"
          write(outunit,*)"STOP file found. Checkpointing and Stopping"
          write(outunit,*)"REMEMBER to remove file before restarting"
          call check(ncycle,time,nvar,naux)
          stop
       endif

      go to 20
c
999   continue

c
c  # computation is complete to final time or requested number of steps
c
       if (ncycle .ge. nstop .and. tfinal .lt. rinfinity) then
c         # warn the user that calculation finished prematurely
          write(outunit,102) nstop
          write(6,102) nstop
  102     format('*** Computation halted after nv(1) = ',i8,
     &           '  steps on coarse grid')
          endif
c
c  # final output (unless we just did it above)
c
      dump_final = ((iout.lt.iinfinity) .and. (mod(ncycle,iout).ne.0))
      if (.not. dumpout) then
          if (nout > 0) then
              dump_final = (tout(nout).eq.tfinal)
              endif
          endif
      
      if (dump_final) then
           call valout(1,lfine,time,nvar,naux)
           if (printout) call outtre(mstart,.true.,nvar,naux)
      endif

c ## tick timing moved here so can be saved in checkpoint file 
c
      call system_clock(tick_clock_finish,tick_clock_rate)
      call cpu_time(tick_cpu_finish)
      timeTick = timeTick + tick_clock_finish - tick_clock_start 
      timeTickCPU = timeTickCPU + tick_cpu_finish - tick_cpu_start 
 

c  # checkpoint everything for possible future restart
c  # (unless we just did it based on dumpchk)
c
      if (checkpt_style .ne. 0) then  ! want a chckpt
         ! check if just did it so dont do it twice
         if (.not. dumpchk) call check(ncycle,time,nvar,naux)
      else  ! no chkpt wanted, so need to print gauges separately
         if (num_gauges .gt. 0) then
            do ii = 1, num_gauges
               call print_gauges_and_reset_nextLoc(ii)
            end do
         endif
      endif
             
   

      write(6,*) "Done integrating to time ",time
      return
      end
