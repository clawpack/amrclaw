c
c ---------------------------------------------------------
c
      subroutine check(nsteps,time,nvar,naux)
c
c :::::::::::::::::::::: CHECK ::::::::::::::::::::::::::::::::;
c   check point routine - can only call at end of coarse grid cycle
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      use amr_module
      use gauges_module, only:  num_gauges
      use gauges_module, only: print_gauges_and_reset_nextLoc

      implicit double precision (a-h,o-z)
      integer tchkunit
      parameter (tchkunit = 13)
      character  chkname*13
      character  tchkname*13

      write(6,601) time,nsteps
 601  format('Creating checkpoint file at t = ',e16.9,'  nsteps = ',i5)
c
      if (checkpt_style < 0) then

c         #  Alternate between two sets of files, overwriting the oldest
c         #  one, so that files do not accumulate with frequent checkpoints.
c
c         # Note that logical check_a is stored in amr_module, initialized
c         # in amr2 and perhaps reset properly in restrt.

          if (check_a) then
              chkname = 'fort.chkaaaaa'
              tchkname = 'fort.tckaaaaa'
            else
              chkname = 'fort.chkbbbbb'
              tchkname = 'fort.tckbbbbb'
            endif
          check_a = .not. check_a   ! to use other file next time

      else

c         # make a new checkpoint file for this time
c         ###  make the file name showing the time step
c
          chkname = 'fort.chkxxxxx'
          tchkname = 'fort.tckxxxxx'
          nstp = nsteps
          do 20 ipos = 13, 9, -1
             idigit = mod(nstp,10)
             chkname(ipos:ipos) = char(ichar('0') + idigit)
             tchkname(ipos:ipos) = char(ichar('0') + idigit)
             nstp = nstp / 10
 20       continue
      endif

      open(unit=tchkunit,file=tchkname,status='unknown',
     .     form='formatted')
      open(unit=chkunit,file=chkname,status='unknown',
     .     form='unformatted')
c
c     ###  dump the data
c
      write(chkunit) lenmax,lendim,memsize
      write(chkunit) (alloc(i),i=1,lendim)
      write(chkunit) hxposs,hyposs,hzposs,possk,icheck
      write(chkunit) lfree,lenf
      write(chkunit) rnode,node,lstart,newstl,listsp,tol,
     1          ibuff,mstart,ndfree,ndfree_bnd,lfine,iorder,mxnest,
     2          intratx,intraty,intratz,kratio,iregsz,jregsz,kregsz,
     3          iregst,jregst,kregst,iregend,jregend,kregend,
     4          numgrids,kcheck,nsteps,time,matlabu
      write(chkunit) avenumgrids, iregridcount,
     1               evol, rvol, rvoll, lentot, tmass0,cflmax
c
      close(chkunit)

c     # flush open running output files fort.amr, fort.gauge, fort.debug
c     # if code dies it will at least have output up to this checkpoint time

      flush(outunit)        ! defined in amr_module.f90
      flush(dbugunit)       ! defined in amr_module.f90
c      flush(OUTGAUGEUNIT)   ! defined in gauges_module.f90
      do ii = 1, num_gauges
         call print_gauges_and_reset_nextLoc(ii,nvar)
      end do

c     # write the time stamp file last so it's not updated until data is
c     # all dumped, in case of crash mid-dump.
      write(tchkunit,*) 'Checkpoint file at time t = ',time
      write(tchkunit,*) 'alloc size memsize = ',memsize
      write(tchkunit,*) 'Number of steps taken = ',nsteps
      close(tchkunit)
c
      return
      end
