c :::::::::::::::::::::::::::: RELOAD ::::::::::::::::::::::::::::::::
c read back in the check point files written by subr. check.
c
c Modified version of restrt that only loads old data
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c ---------------------------------------------------------
c
      subroutine reload(adjfile, k)
c
      use amr_reload_module
      use adjoint_module

      implicit none

      integer, intent(in) :: k
      integer :: i
      character(len=*), intent(in) :: adjfile
      logical foundFile, initial

c     ! Checking to see if file exists
      write(6,*) 'Attempting to reload data '
      write(6,*) '  checkpoint file: ',trim(adjfile)
      inquire(file=trim(adjfile),exist=foundFile)
      if (.not. foundFile) then
        write(*,*)" Did not find checkpoint file!"
        stop
      endif
      open(rstunit,file=trim(adjfile),status='old',form='unformatted')
      rewind rstunit

c     ! Starting the bulk of the reading
      read(rstunit) adjoints(k)%lenmax,adjoints(k)%lendim,
     .    adjoints(k)%isize

      allocate(adjoints(k)%alloc(adjoints(k)%isize))

      memsize = adjoints(k)%isize

      read(rstunit) (adjoints(k)%alloc(i),i=1,adjoints(k)%lendim)
      read(rstunit) adjoints(k)%hxposs,adjoints(k)%hyposs,
     .       adjoints(k)%possk,adjoints(k)%icheck
      read(rstunit) adjoints(k)%lfree,adjoints(k)%lenf
      read(rstunit) adjoints(k)%rnode,adjoints(k)%node,
     1       adjoints(k)%lstart,adjoints(k)%newstl,
     1       adjoints(k)%listsp, adjoints(k)%tol,
     1       adjoints(k)%ibuff,adjoints(k)%mstart,
     1       adjoints(k)%ndfree,adjoints(k)%ndfree_bnd,
     1       adjoints(k)%lfine,
     1       adjoints(k)%iorder,adjoints(k)%mxnest,
     2       adjoints(k)%intratx,adjoints(k)%intraty,
     2       adjoints(k)%kratio,adjoints(k)%iregsz,
     2       adjoints(k)%jregsz, adjoints(k)%iregst,
     2       adjoints(k)%jregst,adjoints(k)%iregend,
     2       adjoints(k)%jregend, adjoints(k)%numgrids,
     3       adjoints(k)%kcheck,adjoints(k)%nsteps,
     3       adjoints(k)%time, adjoints(k)%matlabu
      read(rstunit) adjoints(k)%avenumgrids,
     1              adjoints(k)%iregridcount, adjoints(k)%evol,
     1              adjoints(k)%rvol,adjoints(k)%rvoll,
     1              adjoints(k)%lentot,adjoints(k)%tmass0,
     1              adjoints(k)%cflmax

      close(rstunit)

      write(outunit,100) adjoints(k)%nsteps,adjoints(k)%time
      write(6,100) adjoints(k)%nsteps,adjoints(k)%time
 100  format(/,' Data comes from calculating over ',i5,' steps',
     1        /,'  (time = ',e15.7,')')

      return
      end
