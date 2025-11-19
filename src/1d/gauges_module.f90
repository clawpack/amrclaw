! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::: Parameters, variables, subroutines related to gauges
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Contains:
!   subroutine set_gauges
!     Called initially to read from gauges.data
!   subroutine setbestsrc
!     Called each time regridding is done to determine which patch to 
!     use for interpolating to each gauge location.
!   subroutine update_gauges
!   subroutine print_gauges_and_reset_nextLoc
!
! Note: by default all components of q are printed at each gauge.
!       ASCII precision/format is controlled by FORMAT 100 in
!       print_gauges_and_reset_nextLoc.
!
! Updates for binary gauge output:
!   - New per-gauge file format flag read from gauges.data
!     (1=ascii, 2=binary32, 3=binary64)
!   - Binary files are written as unformatted stream files:
!       record = [level :: INTEGER], [time :: REAL*4 or REAL*8],
!                [q(1:nvar) :: REAL*4 or REAL*8]
!   - A one-line ASCII header is written once when the .bin file
!     is created: "# file_format: binary32" or "binary64"

module gauges_module

    implicit none
    save

    ! 1=ascii, 2=binary32, 3=binary64
    integer, save, allocatable :: gauge_file_format(:)

    logical, private :: module_setup = .false.

    integer, parameter :: OUTGAUGEUNIT=89
    integer :: num_gauges, inum
    real(kind=8), allocatable :: xgauge(:), t1gauge(:), t2gauge(:)
    integer, allocatable, dimension(:) :: mbestsrc, mbestorder, &
         igauge, mbestg1, mbestg2, nextLoc

    integer, parameter :: MAXDATA=1000    !was 1000
    real(kind=8), pointer :: gaugeArray(:,:,:)
    integer, pointer :: levelArray(:,:)

contains

    subroutine set_gauges(restart, nvar, fname)
        use amr_module
        implicit none

        ! Input
        character(len=*), intent(in), optional :: fname
        logical, intent(in) :: restart
        integer, intent(in) :: nvar

        ! Locals
        character(len=256) :: line
        integer :: ios
        integer, parameter :: iunit = 7
        integer :: i, ipos, idigit
        character(len=24) :: fileName
        character(len=24) :: binName
        character(len=15) :: numstr

        if (.not. module_setup) then

            ! Open gauges.data
            if (present(fname)) then
                call opendatafile(iunit,fname)
            else
                call opendatafile(iunit,'gauges.data')
            endif

            read(iunit,*) num_gauges

            allocate(xgauge(num_gauges))
            allocate(t1gauge(num_gauges), t2gauge(num_gauges))
            allocate(mbestsrc(num_gauges), mbestorder(num_gauges))
            allocate(igauge(num_gauges))
            allocate(mbestg1(maxgr), mbestg2(maxgr))

            allocate(nextLoc(num_gauges))
            allocate(gaugeArray(nvar+1,MAXDATA,num_gauges))  ! +1 for time
            allocate(levelArray(MAXDATA,num_gauges))

            do i=1,num_gauges
                read(iunit,*) igauge(i), xgauge(i), t1gauge(i), t2gauge(i)
            enddo



            if (.not. allocated(gauge_file_format)) allocate(gauge_file_format(num_gauges))
            gauge_file_format = 1

            ! Skip comment/blank lines until we hit the numeric formats line
            do
                  read(iunit,'(A)', iostat=ios) line
                  if (ios /= 0) then
                        exit  ! no more lines; keep default (all ASCII)
                  end if
                  if (len_trim(line) == 0) cycle
                  if (line(1:1) == '#' .or. line(1:1) == '!') cycle
                  ! Found a data line: parse num_gauges integers from it
                  read(line, *, iostat=ios) (gauge_file_format(i), i=1,num_gauges)
                  exit
            end do

            close(iunit)


            ! initialize
            mbestsrc = 0
            nextLoc  = 1  ! next location to be filled with gauge info

            ! Create files / headers per gauge
            do i = 1, num_gauges
               inum = igauge(i)
               write (numstr,'(I0.5)') inum
               fileName = 'gauge'//trim(numstr)//'.txt'
               binName  = 'gauge'//trim(numstr)//'.bin'

               select case (gauge_file_format(i))
               case (1)   ! ASCII (existing)
                  if (restart) then
                     open(unit=OUTGAUGEUNIT, file=fileName, status='old',        &
                          position='append', form='formatted')
                  else
                     open(unit=OUTGAUGEUNIT, file=fileName, status='unknown',    &
                          position='append', form='formatted')
                     rewind OUTGAUGEUNIT
                     write(OUTGAUGEUNIT,100) igauge(i), xgauge(i), nvar
100                  format("# gauge_id= ",i0," location=( ",1e15.7," ) num_eqn= ",i2)
                     write(OUTGAUGEUNIT,101)
101                  format("# Columns: level time q(1 ... num_eqn)")
                  endif
                  close(OUTGAUGEUNIT)

               case (2)   ! binary32
                  if (.not. restart) then
                     open(unit=OUTGAUGEUNIT, file=binName, status='unknown', form='formatted')
                     write(OUTGAUGEUNIT,'(A)') '# file_format: binary32'
                     close(OUTGAUGEUNIT)
                  endif

               case (3)   ! binary64
                  if (.not. restart) then
                     open(unit=OUTGAUGEUNIT, file=binName, status='unknown', form='formatted')
                     write(OUTGAUGEUNIT,'(A)') '# file_format: binary64'
                     close(OUTGAUGEUNIT)
                  endif

               end select
            end do

            module_setup = .true.
        end if

    end subroutine set_gauges

! --------------------------------------------------------------------

    subroutine setbestsrc()
!     Called every time grids change, to set the best source grid patch
!     for each gauge, i.e. the finest level patch that includes the gauge.

      use amr_module
      implicit none

      integer :: lev, mptr, i, k1, ki

      ! initialize sources to 0 then check that all set
      do i = 1, num_gauges
         mbestsrc(i) = 0
      end do

      do 20 lev = 1, lfine  
          mptr = lstart(lev)
 5        do 10 i = 1, num_gauges
            if ((xgauge(i) .ge. rnode(cornxlo,mptr)) .and. &
                (xgauge(i) .le. rnode(cornxhi,mptr)) ) then
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

      ! Sort source arrays for easy testing during integration
      call qsorti(mbestorder,num_gauges,mbestsrc)

      ! Figure out the set of gauges that should be handled on each grid:
      mbestg1 = 0
      mbestg2 = 0

      k1 = 0
      do i=1,num_gauges
          ki = mbestsrc(mbestorder(i))
          if (ki > k1) then
              if (k1 > 0) then
                  mbestg2(k1) = i-1
              endif
              mbestg1(ki) = i
          endif
          k1 = ki
      enddo
      if (num_gauges > 0) then
          mbestg2(ki) = num_gauges
      endif

    end subroutine setbestsrc

! -------------------------------------------------------------------------

    subroutine update_gauges(q,aux,xlow,nvar,mitot,naux,mptr)
!
!     Called each time step for each grid patch, to output gauge values
!     for all gauges for which this patch is the best one to use.

      use amr_module
      implicit none

      real(kind=8), intent(in) ::  xlow
      integer, intent(in) ::  nvar,mitot,naux,mptr
      real(kind=8), intent(in) ::  q(nvar,mitot)
      real(kind=8), intent(in) ::  aux(naux,mitot)

      ! locals:
      real(kind=8) :: var(maxvar)
      real(kind=8) :: xcent,xoff,tgrid,hx
      integer :: level,i,iindex,ivar, ii,i1,i2
      integer :: nindex

      if (num_gauges == 0) return

      i1 = mbestg1(mptr)
      i2 = mbestg2(mptr)
      if (i1 == 0) return   ! no gauges for this grid

      tgrid = rnode(timemult,mptr)
      level = node(nestlevel,mptr)
      hx    =  hxposs(level)

      do 10 i = i1,i2
        ii = mbestorder(i)
        if (mptr .ne. mbestsrc(ii)) then
            write(6,*) '*** should not happen... i, ii, mbestsrc(ii), mptr:'
            write(6,*) i, ii, mbestsrc(ii), mptr
            stop
        endif
        if (tgrid.lt.t1gauge(ii) .or. tgrid.gt.t2gauge(ii)) go to 10

        ! linear interpolation to gauge location
        iindex =  int(.5d0 + (xgauge(ii)-xlow)/hx)
        if (iindex .lt. nghost .or. iindex .gt. mitot-nghost) &
          write(*,*)"ERROR in output of Gauge Data "
        xcent  = xlow + (iindex-.5d0)*hx
        xoff   = (xgauge(ii)-xcent)/hx

        do ivar = 1, nvar
           var(ivar) = (1.d0-xoff)*q(ivar,iindex) + xoff*q(ivar,iindex+1)
           if (abs(var(ivar)) .lt. 1.d-90) var(ivar) = 0.d0
        end do

        ! save info for this time
        nindex = nextLoc(ii)
        levelArray(nindex,ii) = level
        gaugeArray(1,nindex,ii) = tgrid
        do ivar = 1, nvar
           gaugeArray(1+ivar,nindex,ii) = var(ivar)
        end do

        nextLoc(ii) = nextLoc(ii) + 1
        if (nextLoc(ii) .gt. MAXDATA) then
          call print_gauges_and_reset_nextLoc(ii,nvar)
        endif
 10   continue

    end subroutine update_gauges

! -------------------------------------------------------------------------

    subroutine print_gauges_and_reset_nextLoc(gaugeNum,nvar)
!
!    Array of gauge data for this gauge reached max capacity (or checkpoint):
!    write buffered data to file, then reset buffer.

      implicit none
      integer :: gaugeNum,nvar,j,inum,k,idigit,ipos,myunit
      integer :: omp_get_thread_num, mythread
      character(len=24) :: fileName
      character(len=24) :: binName
      character(len=15) :: numstr

      ! binary32 temporaries
      integer :: lvl
      real(kind=4) :: t32
      real(kind=4), allocatable :: q32(:)

      ! Build names and unit
      inum = igauge(gaugeNum)
      write (numstr,'(I0.5)') inum
      fileName = 'gauge'//trim(numstr)//'.txt'
      binName  = 'gauge'//trim(numstr)//'.bin'

      mythread = 0
!$    mythread = omp_get_thread_num()
      myunit = OUTGAUGEUNIT + mythread

      select case (gauge_file_format(gaugeNum))

      case (1)   ! ASCII (existing behavior)
         open(unit=myunit, file=fileName, status='old', position='append', form='formatted')
         do j = 1, nextLoc(gaugeNum)-1
            write(myunit,100) levelArray(j,gaugeNum),      &
                              (gaugeArray(k,j,gaugeNum),k=1,nvar+1) ! includes time
         end do
         close(myunit)

      case (2)   ! binary32
         open(unit=myunit, file=binName, status='unknown', form='unformatted', &
              access='stream', position='append')
         allocate(q32(nvar))
         do j = 1, nextLoc(gaugeNum)-1
            lvl = levelArray(j,gaugeNum)
            t32 = real(gaugeArray(1,j,gaugeNum), kind=4)               ! time
            q32 = real(gaugeArray(2:1+nvar,j,gaugeNum), kind=4)        ! q(1:nvar)
            write(myunit) lvl, t32, q32
         end do
         deallocate(q32)
         close(myunit)

      case (3)   ! binary64
         open(unit=myunit, file=binName, status='unknown', form='unformatted', &
              access='stream', position='append')
         do j = 1, nextLoc(gaugeNum)-1
            write(myunit) levelArray(j,gaugeNum),  &
                          gaugeArray(1,j,gaugeNum),        &  ! time (kind=8)
                          gaugeArray(2:1+nvar,j,gaugeNum)     ! q(1:nvar) (kind=8)
         end do
         close(myunit)

      end select

      ! Reset buffer; ASCII format statement lives here
      nextLoc(gaugeNum) = 1
100   format(i5.2, 10e15.7)

    end subroutine print_gauges_and_reset_nextLoc

      ! -------------------------------------------------------------------------
      !
            subroutine flush_all_gauges(nvar)
      !     Write any buffered gauge data that hasn't filled MAXDATA yet.
            implicit none
            integer, intent(in) :: nvar
            integer :: ii

            if (num_gauges <= 0) return

            do ii = 1, num_gauges
            if (nextLoc(ii) > 1) then
                  call print_gauges_and_reset_nextLoc(ii, nvar)
            end if
            end do

            end subroutine flush_all_gauges


end module gauges_module
