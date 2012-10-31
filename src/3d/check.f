c
c ---------------------------------------------------------
c
      subroutine check(nsteps,time,nvar,naux)
c
c :::::::::::::::::::::: CHECK ::::::::::::::::::::::::::::::::;
c   check point routine - can only call at end of coarse grid cycle
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      implicit double precision (a-h,o-z)
      character  chkname*12
      include  "call.i"
c
c     ###  make the file name showing the time step
c
      chkname = 'fort.chkxxxx'
      nstp = nsteps
      do 20 ipos = 12, 9, -1
       idigit = mod(nstp,10)
       chkname(ipos:ipos) = char(ichar('0') + idigit)
       nstp = nstp / 10
 20   continue
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
c
      return
      end
