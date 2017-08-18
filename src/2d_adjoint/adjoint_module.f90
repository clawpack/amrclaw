! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::     Module to define and work with adjoint type
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module adjoint_module
    use amr_module

    type adjointData_type

        real(kind=8), allocatable, dimension(:) :: alloc

        real(kind=8) hxposs(maxlv),hyposs(maxlv)

        integer meqn, ngrids, naux, ndim, nghost, lfine
        real(kind=8) time, xlowvals(maxgr), ylowvals(maxgr)
        integer ncellsx(maxgr), ncellsy(maxgr), loc(maxgr)
        integer gridlevel(maxgr), gridpointer(maxgr)

    end type adjointData_type

    type(adjointData_type), allocatable :: adjoints(:)
    integer :: totnum_adjoints, &
               counter, innerprod_index
    real(kind=8) :: trange_start, trange_final
    character(len=365), allocatable :: adj_files(:)
    logical :: adjoint_flagging

contains

    subroutine read_adjoint_data()

        use amr_reload_module
        implicit none

        ! Function Arguments
        character(len=*), parameter :: adjointfile = 'adjoint.data'
        character(len=200) :: adjoint_output
        logical :: fileExists
        integer :: iunit, k, r
        integer :: fileStatus = 0
        real(kind=8) :: finalT
        real(kind=8) :: t1,t2

        ! Read adjoint specific information

        adjoint_flagging = .true.

        inquire(file=adjointfile, exist=fileExists)
        if (fileExists) then

            iunit = 16
            call opendatafile(iunit,adjointfile)
            read(iunit,*) adjoint_output

            ! time period of interest:
            read(iunit,*) t1
            read(iunit,*) t2

            call set_time_window(t1, t2)

            read(iunit,*) totnum_adjoints
            read(iunit,*) innerprod_index
            allocate(adj_files(totnum_adjoints))

            do 20 k = 1, totnum_adjoints
                read(iunit,*) adj_files(totnum_adjoints + 1 - k)
            20  continue
            close(iunit)

            if (size(adj_files) <= 0) then
                print *, 'Error: no adjoint output files found.'
                stop
            endif

            ! Allocate space for the number of needed checkpoint files
            allocate(adjoints(totnum_adjoints))

            do 50 k = 1, totnum_adjoints
                ! Load checkpoint files
                call reload(adj_files(k),k)

            50 continue

        else
            print *, 'Error, adjoint.data file does not exist.'
        endif

    end subroutine read_adjoint_data

    subroutine set_time_window(t1, t2)

        implicit none

        ! Function Arguments
        real(kind=8), intent(in) :: t1, t2

        trange_final = t2
        trange_start = t1
        write(*,*) "Time range: ", t1, t2

    end subroutine set_time_window


end module adjoint_module
