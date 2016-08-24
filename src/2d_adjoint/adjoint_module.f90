! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::     Module to define and work with adjoint type
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module adjoint_module
    use amr_reload_module

    type adjointData_type

        integer lenmax, lendim, isize, lentot
        real(kind=8), allocatable, dimension(:) :: alloc

        real(kind=8) hxposs(maxlv), hyposs(maxlv),possk(maxlv)
        integer icheck(maxlv)

        ! for space management of alloc array
        integer lfree(lfdim,2),lenf

        real(kind=8) rnode(rsize, maxgr)
        integer node(nsize, maxgr), lstart(maxlv),newstl(maxlv), &
                listsp(maxlv)
        integer ibuff, mstart, ndfree, ndfree_bnd, lfine, iorder, mxnest, &
                intratx(maxlv),intraty(maxlv), kratio(maxlv), &
                iregsz(maxlv),jregsz(maxlv), iregst(maxlv),jregst(maxlv), &
                iregend(maxlv),jregend(maxlv), numgrids(maxlv), kcheck, &
                nsteps, matlabu, iregridcount(maxlv)
        real(kind=8) time, tol, rvoll(maxlv),evol,rvol, &
                     cflmax, avenumgrids(maxlv)

        ! variable for conservation checking
        real(kind=8) tmass0

    end type adjointData_type

    type(adjointData_type), allocatable :: adjoints(:)
    integer :: totnum_adjoints, nvar, naux, &
               counter, innerprod_index
    real(kind=8) :: trange_start, trange_final
    real(kind=8), allocatable :: adj_times(:)
    character(len=365), allocatable :: adj_files(:)
    logical :: adjoint_flagging

contains

    subroutine read_adjoint_data(adjointFolder)

        use amr_reload_module
        implicit none

        ! Function Arguments
        character(len=*), parameter :: adjointfile = 'adjoint.data'
        character(len=*), intent(in) :: adjointFolder
        logical :: fileExists
        integer :: iunit, k, r
        integer :: fileStatus = 0
        real(kind=8) :: finalT

        adjoint_flagging = .true.

        inquire(file=adjointfile, exist=fileExists)
        if (fileExists) then

            ! Reload adjoint files
            call amr2_reload(adjointFolder)
            iunit = 16
            call opendatafile(iunit,adjointfile)

            read(iunit,*) totnum_adjoints
            read(iunit,*) innerprod_index
            allocate(adj_times(totnum_adjoints))
            allocate(adj_files(totnum_adjoints))

            do 20 k = 1, totnum_adjoints
                read(iunit,*) adj_files(totnum_adjoints + 1 - k)
                read(iunit,*) adj_times(totnum_adjoints + 1 - k)
            20  continue
            close(iunit)

            finalT = adj_times(1)

            do 60 k = 1, totnum_adjoints
                ! Reverse times
                adj_times(k) = finalT - adj_times(k)
            60  continue

            ! Allocate space for the number of needed checkpoint files
            allocate(adjoints(totnum_adjoints))

            do 50 k = 1, totnum_adjoints
                ! Load checkpoint files
                call reload(adj_files(k),k)

                ! Reverse times
                adjoints(k)%time = adj_times(k)

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

    end subroutine set_time_window

end module adjoint_module
