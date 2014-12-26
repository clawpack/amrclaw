! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::: Parameters, variables, subroutines related to gauges
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Contains:
!   subroutine set_gauges
!     Called initially to read from gauges.data
!   subroutine setbestsrc
!     Called each time regridding is done to determine which patch to 
!     use for interpolating to each gauge location.
!   subroutine print_gauges
!     Called each time step for each grid patch.
!     Refactored dumpgauge routine to interpolate for all gauges on patch.
!
!     Note: by default all components of q are printed at each gauge.
!     To print something different or a different precision, modify 
!     format statement 100 and/or the write statement that uses it.
!   
! Note: Updated for Clawpack 5.3.0:
!   - the dumpgauge and setbestsrc subroutines have been moved to this module 
!     and the dumpgauge subroutine has been refactored and renamed print_gauges.
!   - dumpgauge.f must be removed from Makefiles.
!   - setbestsrc no longer uses quicksort to sort gauge numbers and
!     print_gauges no longer uses binary search to locate first gauge handled
!     by a grid.  Instead loop over all gauges and skip those not needed.
!     This may cause the gauges to be printed in a different order at some
!     time steps, but values printed should be the same as in previous code.

module gauges_module

    implicit none
    save

    integer, parameter :: OUTGAUGEUNIT=89
    integer :: num_gauges
    real(kind=8), allocatable :: xgauge(:), ygauge(:), t1gauge(:), t2gauge(:)
    integer, allocatable ::  mbestsrc(:), igauge(:)

contains

    subroutine set_gauges(fname)

        use amr_module

        implicit none

        ! Input
        character(len=*), intent(in), optional :: fname

        ! Locals
        integer :: i
        integer, parameter :: iunit = 7

        ! Open file
        if (present(fname)) then
            call opendatafile(iunit,fname)
        else
            call opendatafile(iunit,'gauges.data')
        endif

        read(iunit,*) num_gauges

        allocate(xgauge(num_gauges), ygauge(num_gauges))
        allocate(t1gauge(num_gauges), t2gauge(num_gauges))
        allocate(mbestsrc(num_gauges))
        allocate(igauge(num_gauges))
        
        do i=1,num_gauges
            read(iunit,*) igauge(i),xgauge(i),ygauge(i),t1gauge(i),t2gauge(i)
        enddo

        close(iunit)
        
        ! initialize for starters
        mbestsrc = 0

        ! open file for output of gauge data 
        ! ascii file with format determined by the write(OUTGAUGEUNIT,100)
        ! statement in print_gauges
        open(unit=OUTGAUGEUNIT, file='fort.gauge', status='unknown', &
                                form='formatted')

    end subroutine set_gauges


!
! --------------------------------------------------------------------
!
      subroutine setbestsrc()
!
!     Called every time grids change, to set the best source grid patch
!     for each gauge, i.e. the finest level patch that includes the gauge.
!
!     lbase is grid level that didn't change, but since fine
!     grid may have disappeared, we still have to look starting
!     at coarsest level 1.
!
      use amr_module
      implicit none

      integer :: lev, mptr, i

!
! ##  set source grid for each loc from coarsest level to finest.
! ##  that way finest src grid left and old ones overwritten
! ##  this code uses fact that grids do not overlap

! # for debugging, initialize sources to 0 then check that all set
      do i = 1, num_gauges
         mbestsrc(i) = 0
      end do

 
      do 20 lev = 1, lfine  
          mptr = lstart(lev)
 5        do 10 i = 1, num_gauges
            if ((xgauge(i) .ge. rnode(cornxlo,mptr)) .and. &
                (xgauge(i) .le. rnode(cornxhi,mptr)) .and. &  
                (ygauge(i) .ge. rnode(cornylo,mptr)) .and. &
                (ygauge(i) .le. rnode(cornyhi,mptr)) ) then
               mbestsrc(i) = mptr
            endif
 10       continue

          mptr = node(levelptr, mptr)
          if (mptr .ne. 0) go to 5
 20   continue


      do i = 1, num_gauges
        if (mbestsrc(i) .eq. 0) &
            write(6,*)"ERROR in setting grid src for gauge data",i
      end do

!     Removed this call since easy just to loop over all gauges:
!       sort the source arrays for easy testing during integration
!       call qsorti(mbestorder,num_gauges,mbestsrc)

      end subroutine setbestsrc

!
! -------------------------------------------------------------------------
!
      subroutine print_gauges(q,aux,xlow,ylow,nvar,mitot,mjtot,naux,mptr)
!
!     This routine is called each time step for each grid patch, to output
!     gauge values for all gauges for which this patch is the best one to 
!     use (i.e. at the finest refinement level).  

!     It is called after ghost cells have been filled from adjacent grids
!     at the same level, so bilinear interpolation can be used to 
!     to compute values at any gauge location that is covered by this grid.  

!     The grid patch is designated by mptr.
!     We only want to set gauges i for which mbestsrc(i) == mptr.
!     The array mbestsrc is reset after each regridding to indicate which
!     grid patch is best to use for each gauge.

!     This is a refactoring of dumpgauge.f from Clawpack 5.2 
!     that does not assume list is sorted.  Loops over all gauges rather
!     than doing a binary search.

      use amr_module

      implicit none

      real(kind=8), intent(in) ::  q(nvar,mitot,mjtot)
      real(kind=8), intent(in) ::  aux(naux,mitot,mjtot)
      real(kind=8), intent(in) ::  xlow,ylow
      integer, intent(in) ::  nvar,mitot,mjtot,naux,mptr

      ! local variables:
      real(kind=8) :: var(maxvar)
      real(kind=8) :: xcent,ycent,xoff,yoff,tgrid,hx,hy
      integer :: level,i,j,ioff,joff,iindex,jindex,ivar

!     write(*,*) 'in print_gauges with num_gauges, mptr = ',num_gauges,mptr
!     write(*,*) 'mbestscr = ',(mbestsrc(j),j=1,num_gauges)

      if (num_gauges == 0) then
         return
      endif


!     # this stuff the same for all gauges on this grid
      tgrid = rnode(timemult,mptr)
      level = node(nestlevel,mptr)
      hx    =  hxposs(level)
      hy    =  hyposs(level)

!     write(*,*) 'tgrid = ',tgrid

      do 10 i = 1, num_gauges
        if (mptr .ne. mbestsrc(i)) go to 10  ! this patch not used
        if (tgrid.lt.t1gauge(i) .or. tgrid.gt.t2gauge(i)) then
!          # don't output at this time for gauge i
           go to 10
           endif
!
!    ## if we did not skip to line 10, we need to output gauge i:
!    ## prepare to do bilinear interp at gauge location to get vars
!
!    *** Note: changed 0.5 to  0.5d0 etc. ****************************
!
        iindex =  int(.5d0 + (xgauge(i)-xlow)/hx)
        jindex =  int(.5d0 + (ygauge(i)-ylow)/hy)
        if ((iindex .lt. nghost .or. iindex .gt. mitot-nghost) .or. &
            (jindex .lt. nghost .or. jindex .gt. mjtot-nghost)) &
          write(*,*)"ERROR in output of Gauge Data "
        xcent  = xlow + (iindex-.5d0)*hx
        ycent  = ylow + (jindex-.5d0)*hy
        xoff   = (xgauge(i)-xcent)/hx
        yoff   = (ygauge(i)-ycent)/hy
        if (xoff .lt. 0.d0 .or. xoff .gt. 1.d0 .or. &
            yoff .lt. 0.d0 .or. yoff .gt. 1.d0) then
           write(6,*)" BIG PROBLEM in DUMPGAUGE", i
        endif

!       ## bilinear interpolation
        do ivar = 1, nvar
           var(ivar) = (1.d0-xoff)*(1.d0-yoff)*q(ivar,iindex,jindex) &
                   + xoff*(1.d0-yoff)*q(ivar,iindex+1,jindex) &
                   + (1.d0-xoff)*yoff*q(ivar,iindex,jindex+1) &
                   + xoff*yoff*q(ivar,iindex+1,jindex+1)
!          # for printing without underflow of exponent:
           if (abs(var(ivar)) .lt. 1.d-90) var(ivar) = 0.d0
        end do


!$OMP CRITICAL (gaugeio)

!       # output values at gauge, along with gauge no, level, time:
!       # if you want to print out something different at each gauge,
!       # modify this...
        write(OUTGAUGEUNIT,100)igauge(i),level, tgrid,(var(j),j=1,nvar)

!       # if you want to modify number of digits printed, modify this...
100     format(2i5,15e15.7)

!$OMP END CRITICAL (gaugeio)


 10     continue  ! end of loop over all gauges
 
      end subroutine print_gauges

end module gauges_module
