c
c ----------------------------------------------------------------
c
       program amr3ez
c
c  Use adaptive mesh refinement to solve the hyperbolic 3-d equation:
c
c              u  +  f(u)    + g(u)   + h(u)   = 0
c               t         x        y        z
c
c or the more general non-conservation law form:
c              u  +  A u     + B u    + C u    = 0
c               t         x        y       z
c
c  using the wave propagation method as in CLAWPACK in combination
c  with the locally uniform embedded grids of AMR.

c  Estimate error with Richardson extrap. (in errest.f)
c  + gradient checking (in errsp.f).  Initial conditions set
c  in (qinit.f), b.c.'s in (bc3amr.f).

c  Specify rectangular domain from
c           (xlower,ylower,zlower) to (xupper,yupper,zupper)
c
c  No rotated rectangles are used in this version.
c  Periodic b.c.'s finally implemented.
c
c =========================================================================
c  Copyright 2001,  Marsha J. Berger and Randall J. LeVeque
c
c  This software is made available for research and instructional use only.
c  You may copy and use this software without charge for these non-commercial
c  purposes, provided that the copyright notice and associated text is
c  reproduced on all copies.  For all other uses (including distribution of
c  modified versions), please contact the author at the address given below.
c
c  *** This software is made available "as is" without any assurance that it
c  *** will work for your purposes.  The software may in fact have defects, so
c  *** use the software at your own risk.
c
c  --------------------------------------
c    AMRCLAW Version 0.4,  February, 2001
c     compatible with CLAWPACK Version 4.0
c    Homepage: http://www.amath.washington.edu/~claw/
c  --------------------------------------
c
c   Authors:
c
c             Marsha J. Berger
c             Courant Institute of Mathematical Sciences
c             New York University
c             251 Mercer St.
c             New York, NY 10012
c             berger@cims.nyu.edu
c
c             Randall J. LeVeque
c             Applied Mathematics
c             Box 352420
c             University of Washington,
c             Seattle, WA 98195-2420
c             rjl@amath.washington.edu
c
c =========================================================================
c


c
c ----------------------------------------------------------------
c
      use amr_module
      implicit double precision (a-h,o-z)

      common /combc3/ mthbc(6)

      character * 12     pltfile,infile,outfile,dbugfile,matfile
      character*10       matname2
      logical            vtime,rest
      dimension          tout(maxout)

      integer oldmode,omp_get_max_threads
      integer clock_start, clock_finish, clock_rate
      
c
c
c  you may want to turn this on for SUN workstation, or replace
c  set to signal on overflow, divide by zero, and illegal operation
c
c       oldmode = ieee_handler("set","common",SIGFPE_ABORT)
c       if (oldmode .ne. 0) then
c           write(outunit,*)' could not set ieee trapper '
c           write(*,*)      ' could not set ieee trapper '
c           stop
c        endif
c
      infile   = 'amr3ez.data'
      outfile  = 'fort.amr'
      pltfile  = 'fort.ncar'
      dbugfile = 'fort.debug'
      matfile  = 'fort.nplot'

      open(inunit,  file=infile,status='old',form='formatted')
      open(outunit, file=outfile,status='unknown',form='formatted')
      open(dbugunit,file=dbugfile,status='unknown',form='formatted')
c
c     domain variables
      read(inunit,*) nx
      read(inunit,*) ny
      read(inunit,*) nz
      read(inunit,*) mxnest
      if (abs(mxnest) .gt. maxlv) then
         write(outunit,*)
     &    'Error ***   mxnest > max. allowable levels (maxlv) in common'
         write(*,*)
     &    'Error ***   mxnest > max. allowable levels (maxlv) in common'
         stop
      endif
      if (mxnest .lt. 0) then
c         # mxnext<0 flags anisotropic refinement - read in refinement
c         # ratios in x,y,t from next three lines:
          mxnest = -mxnest
          read(inunit,*) (intratx(i),i=1,max(1,mxnest-1))
          read(inunit,*) (intraty(i),i=1,max(1,mxnest-1))
          read(inunit,*) (intratz(i),i=1,max(1,mxnest-1))
          read(inunit,*) (kratio(i), i=1,max(1,mxnest-1))
      else
          read(inunit,*) (intratx(i),i=1,max(1,mxnest-1))
          do i=1,max(1,mxnest-1)
             intraty(i) = intratx(i)
             intratz(i) = intratx(i)
             kratio(i)  = intratx(i)
          enddo
      endif


      read(inunit,*) nout
      if (nout .gt. maxout) then
         write(6,*) 'Error *** need to increase maxout in call.i'
         stop
      endif
      read(inunit,*) outstyle
      if (outstyle.eq.1) then
         read(inunit,*) tfinal
         iout = 0
c        # array tout is set below after reading t0
      endif
      if (outstyle.eq.2) then
         read(inunit,*) (tout(i), i=1,nout)
         iout = 0
      endif
      if (outstyle.eq.3) then
         read(inunit,*) iout,nstop
         nout = 0
      endif
      read(inunit,*) possk(1)
      read(inunit,*) dt_max
      read(inunit,*) cflv1
      read(inunit,*) cfl
      read(inunit,*) nv1
      if (outstyle.eq.1 .or. outstyle.eq.2) then
         nstop = nv1
      endif

      read(inunit,*) method(1)
      vtime = (method(1) .eq. 1)
      read(inunit,*) method(2)
      iorder = method(2)
      read(inunit,*) method(3)
      if (method(3) .lt. 0) then
         write(6,*) '*** ERROR ***  method(3) < 0'
         write(6,*) '    dimensional splitting not supported in amrclaw'
         stop
      endif

      read(inunit,*) method(4)
      read(inunit,*) method(5)
      read(inunit,*) mcapa1
      read(inunit,*) naux
      if (naux .gt. maxaux) then
         write(outunit,*) 'Error ***   naux > maxaux in common'
         write(*,*)       'Error ***   naux > maxaux in common'
         stop
      endif
      do iaux = 1, naux
        read(inunit,*) auxtype(iaux)
      end do


      read(inunit,*) nvar
      read(inunit,*) mwaves
      if (mwaves .gt. maxwave) then
         write(outunit,*) 'Error ***   mwaves > maxwave in common'
         write(*,*)       'Error ***   mwaves > maxwave in common'
         stop
      endif
      read(inunit,*) (mthlim(mw), mw=1,mwaves)


      read(inunit,*) t0
      tstart = t0                 ! for common block
      read(inunit,*) xlower
      read(inunit,*) xupper
      read(inunit,*) ylower
      read(inunit,*) yupper
      read(inunit,*) zlower
      read(inunit,*) zupper

      read(inunit,*) nghost
      read(inunit,*) mthbc(1)
      read(inunit,*) mthbc(2)
      read(inunit,*) mthbc(3)
      read(inunit,*) mthbc(4)
      read(inunit,*) mthbc(5)
      read(inunit,*) mthbc(6)

      xperdom = (mthbc(1).eq.2 .and. mthbc(2).eq.2)
      yperdom =  (mthbc(3).eq.2 .and. mthbc(4).eq.2)
      zperdom =  (mthbc(5).eq.2 .and. mthbc(6).eq.2)

      if ((mthbc(1).eq.2 .and. mthbc(2).ne.2) .or.
     &    (mthbc(2).eq.2 .and. mthbc(1).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions: '
         write(6,*) '  mthbc(1) and mthbc(2) must BOTH be set to 2'
         stop
      endif

      if ((mthbc(3).eq.2 .and. mthbc(4).ne.2) .or.
     &    (mthbc(4).eq.2 .and. mthbc(3).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions: '
         write(6,*) '  mthbc(3) and mthbc(4) must BOTH be set to 2'
         stop
      endif

      if ((mthbc(5).eq.2 .and. mthbc(6).ne.2) .or.
     &    (mthbc(6).eq.2 .and. mthbc(5).ne.2)) then
         write(6,*) '*** ERROR ***  periodic boundary conditions: '
         write(6,*) '  mthbc(5) and mthbc(6) must BOTH be set to 2'
         stop
      endif

      if (outstyle.eq.1) then
	   do i=1,nout
	      tout(i) = t0 + i*(tfinal-t0)/float(nout)
	      enddo
           endif

c     restart and checkpointing
      read(inunit,*) rest
      read(inunit,*) ichkpt
c
c     refinement variables
      read(inunit,*) tol
      read(inunit,*) tolsp
      read(inunit,*) kcheck
      read(inunit,*) ibuff
      read(inunit,*) cut
c
c     style of output
c
      read(inunit,*) printout
      read(inunit,*) ncarout
      read(inunit,*) matlabout
c
c
c     # read verbose/debugging flags
      read(inunit,*) dprint
      read(inunit,*) eprint
      read(inunit,*) edebug
      read(inunit,*) gprint
      read(inunit,*) nprint
      read(inunit,*) pprint
      read(inunit,*) rprint
      read(inunit,*) sprint
      read(inunit,*) tprint
      read(inunit,*) uprint
c
c     # close input file
      close(inunit)

c     # look for capacity function via auxtypes:
      mcapa = 0
      do 20 iaux = 1, naux
        if (auxtype(iaux) .eq. "capacity") then
          if (mcapa .ne. 0) then
            write(*,*)" only 1 capacity allowed"
            stop
          else
            mcapa = iaux
          endif
        endif

c       # change to Version 4.1 terminology:
        if (auxtype(iaux) .eq. "leftface") auxtype(iaux) = "xleft"
        if (auxtype(iaux) .eq. "frontface") auxtype(iaux) = "yleft"
        if (auxtype(iaux) .eq. "bottomface") auxtype(iaux) = "zleft"

        if (.not. (auxtype(iaux) .eq."xleft" .or.
     .             auxtype(iaux) .eq. "yleft".or.
     .             auxtype(iaux) .eq. "zleft".or.
     .             auxtype(iaux) .eq. "capacity".or.
     .             auxtype(iaux) .eq. "center"))  then
                  write(*,*)" unknown type for auxiliary variables"
                  write(6,*) auxtype(iaux)
                  write(*,*)" use  xleft/yleft/zleft/center/capacity"
                  stop
        endif
   20   continue

c     ::: error checking of input data :::::::::::::::::::::::

      if (mcapa .ne. mcapa1) then
         write(outunit,*) 'Error ***  mcapa does not agree with auxtype'
         write(*,*) 'Error ***  mcapa does not agree with auxtype'
         stop
         endif


      if (nvar .gt. maxvar) then
         write(outunit,*) 'Error ***   nvar > maxvar in common'
         write(*,*)       'Error ***   nvar > maxvar in common'
         stop
         endif
      if (2*nghost .gt. min(nx,ny,nz)) then
         mindim = 2*nghost
         write(outunit,*) 'Error ***   need finer domain >',
     .         mindim, ' cells'
         write(*,*)       'Error ***   need finer domain >',
     .         mindim, ' cells'
         stop
         endif
      if (mcapa .gt. naux) then
         write(outunit,*) 'Error ***   mcapa > naux in input file'
         write(*,*)       'Error ***   mcapa > naux in input file'
         stop
         endif
      if (.not. vtime .and. nout .ne. 0) then
         write(outunit,*)' cannot specify output times with fixed dt'
         write(*,*)      ' cannot specify output times with fixed dt'
         stop
         endif
c
c
      if (ncarout)
     .   open(pltunit1,file=pltfile,status='unknown',form='formatted')
c

c     # call user routine to set up problem parameters:

      call setprob()
c
      matlabu   = 0
      hxposs(1) = (xupper - xlower) / nx
      hyposs(1) = (yupper - ylower) / ny
      hzposs(1) = (zupper - zlower) / nz
      cflmax = 0.d0
c
c
c
      if (rest) then
         call restrt(nsteps,time,nvar)
         nstart  = nsteps
         write(6,*) ' '
         write(6,*) 'Restarting from previous run'
         write(6,*) '   at time = ',time
         write(6,*) ' '
      else
         lentot = 0
         lenmax = 0
         lendim = 0
         rvol   = 0.0d0
         do 8 i   = 1, mxnest
 8         rvoll(i) = 0.0d0
         evol   = 0.0d0
         call   stst1
         call   domain (nvar,vtime,nx,ny,nz,naux,t0)
         call   setgrd (nvar,cut,naux,dtinit,t0)       
         if (possk(1) .gt. dtinit*cflv1/cfl .and. vtime) then
c        ## initial time step was too large. reset to dt from setgrd
              write(6,*) "*** Initial time step reset for desired cfl"
              possk(1) = dtinit
              do i = 2, mxnest-1
                 possk(i) = possk(i-1)*kratio(i-1)
              end do
         endif

         time = t0
         nstart = 0
      endif

!$    write(outunit,*)" max threads set to ",omp_get_max_threads()
!$    write(*,*)" max threads set to ",omp_get_max_threads()
c
c  print out program parameters for this run
c
      write(outunit,107)tol,tolsp,iorder,kcheck,ibuff,nghost,cut,
     1            mxnest,ichkpt,cfl
      write(outunit,109) xupper,yupper,zupper,
     &      xlower,ylower,zlower,nx,ny, nz
      write(outunit,139)(intratx(i),i=1,mxnest)
      write(outunit,139)(intraty(i),i=1,mxnest)
      write(outunit,139)(intratz(i),i=1,mxnest)
      write(outunit,119) naux
      write(outunit,129) (iaux,auxtype(iaux),iaux=1,naux)
      if (mcapa .gt. 0) write(outunit,149) mcapa
107   format(/
     *       ' amrclaw parameters:',//,
     *       ' error tol            ',e12.5,/,
     *       ' spatial error tol    ',e12.5,/,
     *       ' order of integrator     ',i9,/,
     *       ' error checking interval ',i9,/,
     *       ' buffer zone size        ',i9,/,
     *       ' nghost                  ',i9,/,
     *       ' volume ratio cutoff  ',e12.5,/,
     *       ' max. refinement level   ',i9,/,
     *       ' user sub. calling times ',i9,/,
     *       ' cfl # (if var. delt) ',e12.5,/)
109   format(' xupper(upper corner) ',e12.5,/,
     *       ' yupper(upper corner) ',e12.5,/,
     *       ' zupper(upper corner) ',e12.5,/,
     *       ' xlower(lower corner) ',e12.5,/,
     *       ' ylower(lower corner) ',e12.5,/,
     *       ' zlower(lower corner) ',e12.5,/,
     *       ' nx = no. cells in x dir.',i9,/,
     *       ' ny = no. cells in y dir.',i9,/,
     *       ' nz = no. cells in z dir.',i9,/,/)
139   format(' refinement ratios:      ',6i5,/)
119   format(' no. auxiliary vars.     ',i9)
129   format('       var ',i5,' of type ', a10)
149   format(' capacity fn. is aux. var',i9)
c
      write(6,*) ' '
      write(6,*) 'running amr3ez...  '
      write(6,*) ' '

      call outtre (mstart,printout,nvar,naux)
      write(outunit,*) "  original total mass ..."
      call conck(1,nvar,naux,time,rest)
      call valout(1,lfine,time,nvar,naux)

      ! Timing 
      call system_clock(clock_start,clock_rate)
c
c     --------------------------------------------------------
c     # tick is the main routine which drives the computation:
c     --------------------------------------------------------
      call tick(nvar,iout,nstart,nstop,cut,vtime,time,ichkpt,naux,
     &          nout,tout,rest,dt_max)
c     --------------------------------------------------------

      call system_clock(clock_finish,clock_rate)
      write(outunit,800) dble(clock_finish-clock_start)
     &                  /dble(clock_rate)
 800  format("Total time to solution = ",f16.8," s")

      do level = 1, mxnest
         write(outunit,801) level,tvoll(level)/clock_rate
         write(*,801) level,tvoll(level)/clock_rate
 801     format("Total advanc time on level ",i3," = ",f16.8," s")
      end do

c     # Done with computation, cleanup:

      lentotsave = lentot
      call cleanup(nvar,naux)
      if (lentot .ne. 0) then
        write(outunit,*) lentot," words not accounted for ",
     &                   "in memory cleanup"
        write(*,*)        lentot," words not accounted for ",
     &                   "in memory cleanup"
      endif
c
c report on statistics
c
      open(matunit,file=matfile,status='unknown',form='formatted')
      write(matunit,*) matlabu-1
      write(matunit,*) mxnest
      close(matunit)

      write(outunit,*)
      write(outunit,*)
      do i = 1, mxnest
        if (iregridcount(i) > 0) then
          write(outunit,802) i,avenumgrids(i)/iregridcount(i),
     1                       iregridcount(i)
 802      format("for level ",i3, " average num. grids = ",f10.2,
     1           " over ",i10," steps")
        endif
      end do

      write(outunit,*)
      write(outunit,*)
      write(outunit,901) lentotsave
      write(outunit,902) lenmax
      write(outunit,903) lendim

      write(outunit,904) rvol
      do 60 level = 1,mxnest
 60     write(outunit,905) level, rvoll(level)

      write(outunit,906) evol
      if (evol+rvol .gt. 0.) then
         ratmet = rvol / (evol+rvol) * 100.0d0
      else
         ratmet = 0.0d0
      endif
      write(outunit,907) ratmet
      write(outunit,908) cflmax

 901  format('current  space usage = ',i12)
 902  format('maximum  space usage = ',i12)
 903  format('need space dimension = ',i12,/)
 904  format('number of cells advanced for time integration = ',f20.6)
 905  format(3x,'# cells advanced on level ',i4,' = ',f20.2)
 906  format('number of cells advanced for error estimation = ',f20.6,/)
 907  format(' percentage of cells advanced in time  = ', f10.2)
 908  format(' maximum Courant number seen = ', f10.2)
c
      write(outunit,909)
 909  format(//,' ------  end of AMRCLAW integration --------  ')
c
c     # Close output and debug files.
      close(outunit)
      close(dbugunit)
c
      stop
      end
