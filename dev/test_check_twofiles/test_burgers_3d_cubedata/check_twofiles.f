c
c ---------------------------------------------------------
c
      subroutine check(nsteps,time,nvar,naux)
c
c :::::::::::::::::::::: CHECK ::::::::::::::::::::::::::::::::;
c   check point routine - can only call at end of coarse grid cycle
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      use amr_module
      implicit double precision (a-h,o-z)
      integer tchkunit
      parameter (tchkunit = 13)
      character  chkname*13
      character  tchkname*13
      logical check_a
      common /check_switch/ check_a

      write(6,601) time,nsteps
 601  format('Creating checkpoint file at t = ',e16.9,'  nsteps = ',i5)
c
c     #  Alternate between two sets of files, overwriting the oldest
c     #  one, so that files do not accumulate when doing frequent checkpoints.
c
c     # Note that logical check_a is currently random at start, should
c     # initialize properly when this is rewritten as a module.

      if (check_a) then
          chkname = 'fort.chkaaaaa'
          tchkname = 'fort.tckaaaaa'
        else
          chkname = 'fort.chkbbbbb'
          tchkname = 'fort.tckbbbbb'
        endif
      check_a = .not. check_a   ! to use other file next time

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
     1          ibuff,mstart,ndfree,lfine,iorder,mxnest,
     2          intratx,intraty,intratz,kratio,iregsz,jregsz,kregsz,
     3          iregst,jregst,kregst,iregend,jregend,kregend,
     4          numgrids,kcheck,nsteps,time,matlabu
      write(chkunit) avenumgrids, iregridcount,
     1               evol, rvol, rvoll, lentot, tmass0,cflmax
c
      close(chkunit)

c     # write the time stamp file last so it's not updated until data is
c     # all dumped, in case of crash mid-dump.
      write(tchkunit,*) 'Checkpoint file at time t = ',time
      write(tchkunit,*) 'alloc size memsize = ',memsize
      write(tchkunit,*) 'Number of steps taken = ',nsteps
      close(tchkunit)
c
      return
      end
