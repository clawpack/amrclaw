! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::: Parameters and variables related to gauges
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
module gauges_module

    implicit none
    save

    integer, parameter :: OUTGAUGEUNIT=89
    integer :: num_gauges
    real(kind=8), allocatable :: xgauge(:), ygauge(:), t1gauge(:), t2gauge(:)
    integer, allocatable ::  mbestsrc(:), mbestorder(:), igauge(:)

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
        allocate(mbestsrc(num_gauges), mbestorder(num_gauges))
        allocate(igauge(num_gauges))
        
        do i=1,num_gauges
            read(iunit,*) igauge(i),xgauge(i),ygauge(i),t1gauge(i),t2gauge(i)
        enddo

        close(iunit)
        
        ! initialize for starters
        mbestsrc = 0

        ! open file for output of gauge data all data is output in one binary 
        ! file with format gauge number, level, time, depth by dumpgauge.
        open(unit=OUTGAUGEUNIT, file='fort.gauge', status='unknown', &
                                form='formatted')

    end subroutine set_gauges


end module gauges_module