c=========================================================================
      subroutine setgauges
c=========================================================================

      use amr_module
	  use gauges_module
      implicit double precision (a-h,o-z)
      character*25 fname
      logical foundFile


      fname  = 'gauges.data'
      inquire(file=fname,exist=foundFile)
      if (.not. foundFile) then
        write(*,*) 'You must provide a file ', fname 
        stop
      endif

      iunit = 7
      call opendatafile(iunit,fname)

      read(iunit,*) mgauges
      write(6,*) '+++ in setgauges: mgauges = ',mgauges
      if (mgauges.gt.maxgauges) then
            write(*,*) 'ERROR in setgauges'
            write(*,*) 'mgauges = ',mgauges,'   maxgauges = ',maxgauges
            write(*,*) 'increase maxgauges in amr_module.f90'
            stop
            endif

      do i=1,mgauges
          read(iunit,*) igauge(i),xgauge(i),ygauge(i),
     &                  t1gauge(i),t2gauge(i)
          mbestsrc(i) = 0   ! initialize for starters
          enddo
      close(iunit)

c     # open file for output of gauge data
c     # all data is output in one binary file with format
c     #   gauge number, level, time, depth
c     # by dumpgauge.

      open(unit=OUTGAUGEUNIT,file='fort.gauge',status='unknown',
     .           form='formatted')

      return
      end
