! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::     Module to define and work with adjoint type
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module adjoint_module
    use amr_module

    type adjointData_type

        real(kind=8), allocatable, dimension(:) :: alloc

        real(kind=8) hxposs(maxlv)

        integer meqn, ngrids, naux, ndim, nghost, lfine
        real(kind=8) time, xlowvals(maxgr)
        integer ncellsx(maxgr), loc(maxgr)
        integer gridlevel(maxgr), gridpointer(maxgr)

        ! variable for conservation checking
        real(kind=8) tmass0

    end type adjointData_type

    type(adjointData_type), allocatable :: adjoints(:)
    integer :: totnum_adjoints, &
               counter, innerprod_index
    real(kind=8) :: trange_start, trange_final, levtol(maxlv)
    character(len=365), allocatable :: adj_files(:)
    logical :: adjoint_flagging
    real(kind=8), allocatable, dimension(:) :: errors
    integer, allocatable, dimension(:) :: eptr

contains

    subroutine read_adjoint_data()

        implicit none

        ! Function Arguments
        character(len=*), parameter :: adjointfile = 'adjoint.data'
        character(len=400) :: adjoint_output
        logical :: fileExists
        integer :: iunit, k
        real(kind=8) :: t1,t2

        adjoint_flagging = .true.
        levtol = 0.0d0

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
            print *, 'Error: adjoint.data file does not exist.'
        endif

    end subroutine read_adjoint_data

    subroutine set_time_window(t1, t2)

        implicit none

        ! Function Arguments
        real(kind=8), intent(in) :: t1, t2

        trange_final = t2
        trange_start = t1

    end subroutine set_time_window

    subroutine calculate_tol(lcheck)

       use amr_module
       implicit none

       ! Function Arguments
       integer, intent(in) :: lcheck
       real(kind=8) :: errtotal, cutoff
       real(kind=8) :: sorted(numcells(lcheck)), dt
       integer :: celln

       ! Setting our goal for the maximum amount of error 
       ! for this level
       dt = possk(lcheck)
       cutoff = tol*dt/tfinal ! Total error allowed in this time step
       cutoff = cutoff / 2**(lcheck)

       ! Sorting errors
       call qsortr(sorted, numcells(lcheck), errors)

       errtotal = 0
       do celln = 1, numcells(lcheck)
           errtotal = errtotal + errors(sorted(celln))
           if (errtotal .ge. cutoff) then
               levtol(lcheck) = errors(sorted(celln-1))
               EXIT
           endif
       end do

    end subroutine calculate_tol

end module adjoint_module
