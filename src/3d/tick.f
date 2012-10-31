c
c  -------------------------------------------------------------
c
      subroutine tick(nvar,iout,nstart,nstop,cut,vtime,time,ichkpt,
     &                naux,nout,tout,rest)
c
      implicit double precision (a-h,o-z)

      include  "call.i"

      logical    vtime, dumpout
      dimension dtnew(maxlv), ntogo(maxlv), tlevel(maxlv)
      dimension  tout(nout)
c
c :::::::::::::::::::::::::::: TICK :::::::::::::::::::::::::::::
c  main driver routine.  controls:
c        integration  of all grids.
c        error estimation / regridding
c        output counting
c        updating of fine to coarse grids

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
      if ( iout .eq. 0) iout  = iinfinity
      if (ichkpt .eq. 0) ichkpt = iinfinity
      if (nout .gt. 0) then
	 tfinal = tout(nout)
c        if this is a restart, make sure output times start after restart time
	 if (nstart .gt. 0) then
	    do ii = 1, nout
	      if (tout(ii) .gt. time) then
	        nextout = ii
	        go to 2
	      endif
	    end do
	 else
	    nextout   = 1
	 endif
      else
	 tfinal  = rinfinity
	 nextout = 1     ! keeps it out of the way
      endif
 2    tlevel(1)      = time
      do 5 i       = 2, mxnest
       tlevel(i) = tlevel(1)
 5     continue
c
c  ------ start of coarse grid integration loop. ------------------
c
 20   if (ncycle .ge. nstop .or. time .ge. tfinal) goto 999

      if (nextout  .le. nout) then
         outtime       = tout(nextout)
      else
         outtime       = rinfinity
      endif

      if (time + possk(1) .ge. outtime) then
        if (vtime) then
c          ## adjust time step  to hit outtime exactly, and make output
           possk(1) = outtime - time
           do 12 i = 2, mxnest
 12           possk(i) = possk(i-1) / kratio(i-1)
        endif
c          ## if not variable time step, just turn on output - within tol.
        nextout = nextout + 1
        dumpout = .true.
      else
	dumpout = .false.
      endif
c  take output stuff from here - put it back at end.
c
          level        = 1
	  ntogo(level) = 1
	  do 10 i = 1, maxlv
 10          dtnew(i)  = rinfinity
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
          call regrid(nvar,lbase,cut,naux)
c	  call conck(1,nvar,time,rest)
c         call outtre(lstart(lbase+1),.true.,nvar,naux)
c         call valout(lbase,lfine,tlevel(lbase),nvar,naux)
c
c  maybe finest level in existence has changed. reset counters.
c
          if (rprint .and. lbase .lt. lfine) then
             call outtre(lstart(lbase+1),.false.,nvar,naux)
          endif
 70       continue
          do 80  i  = lbase, lfine
 80          icheck(i) = 0
          do 81  i  = lbase+1, lfine
 81          tlevel(i) = tlevel(lbase)
c
c  ------- done regridding --------------------
c
c integrate all grids at level 'level'.
c
 90       continue

          call advanc(level,nvar,dtlevnew,vtime,naux)

c         # rjl modified 6/17/05 to print out *after* advanc and print cfl
          timenew = tlevel(level)+possk(level)
100       format(' AMRCLAW3... level ',i2,'  CFL = ',e10.4,
     &           '  dt = ',e10.4,  '  t = ',e10.4)
c100       format(' integrating grids at level ',i3,' from t =',
c    &           e11.5,  '  using dt = ',e11.5)

          if (tprint) then
              write(outunit,100)level,cflmax,possk(level),timenew
              endif
          if (method(4).ge.level) then
              write(6,100)level,cflmax,possk(level),timenew
              endif



c        # to debug individual grid updates...
c        call valout(1,lfine,time,nvar,naux)
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
 106		 if ((possk(level)-dtnew(level))/dtnew(level)
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
                 call update(level,nvar)
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
	  call conck(1,nvar,time,rest)

       if (mod(ncycle,ichkpt).eq.0) then
	  call check(ncycle,time,nvar,naux)
       endif
       if ((mod(ncycle,iout).eq.0) .or. dumpout) then
          call valout(1,lfine,time,nvar,naux)
          if (printout) call outtre(mstart,.true.,nvar,naux)
c          write(6,*) ncycle,iout,nvar,time,possk(1)
c          if (visout) call vizout(time,nvar)

       endif

      if ( vtime) then
c
c         find new dt for next cycle (passed back from integration routine).
         do 115 i = 2, lfine
	   ii = lfine+1-i
	         dtnew(ii) = min(dtnew(ii),dtnew(ii+1)*kratio(ii))
 115     continue
         possk(1) = dtnew(1)
         do 120 i = 2, mxnest
 120       possk(i) = possk(i-1) / kratio(i-1)

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
      if (.not. dumpout) then
         if (mod(ncycle,iout).ne.0) then
           call valout(1,lfine,time,nvar,naux)
           if (printout) call outtre(mstart,.true.,nvar,naux)
         endif
      endif

c  # checkpoint everything for possible future restart
c  # (unless we just did it based on ichkpt)
c
c  # don't checkpoint at all if user set ichkpt=0

      if (ichkpt .lt. iinfinity) then
         if ((ncycle/ichkpt)*ichkpt .ne. ncycle) then
           call check(ncycle,time,nvar,naux)
         endif
       endif

      return
      end
