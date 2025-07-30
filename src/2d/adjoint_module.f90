! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::     Module to define and work with adjoint type
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module adjoint_module
    use amr_module

    type adjointData_type

        real(kind=8), allocatable, dimension(:) :: alloc

        real(kind=8) hxposs(maxlv),hyposs(maxlv)

        integer meqn, ngrids, naux, ndim, nghost, lfine
        real(kind=8) time
        real(kind=8),allocatable,dimension(:) :: xlowvals,ylowvals
        integer, allocatable, dimension(:) :: ncellsx,ncellsy,loc, &
                gridlevel,gridpointer

    end type adjointData_type

    type(adjointData_type), allocatable :: adjoints(:)
    integer :: totnum_adjoints, &
               counter, innerprod_index
    real(kind=8) :: trange_start, trange_final, levtol(maxlv)
    character(len=365), allocatable :: adj_files(:)
    logical :: adjoint_flagging
    real(kind=8), allocatable, dimension(:) :: errors
    integer, allocatable, dimension(:) :: eptr
    integer, allocatable, dimension(:) :: grid_num

contains

! ========================================================================
!  read_adjoint_data()
! ========================================================================
    subroutine read_adjoint_data()

        implicit none

        ! Function Arguments
        character(len=*), parameter :: adjointfile = 'adjoint.data'
        character(len=400) :: adjoint_output
        logical :: fileExists
        integer :: iunit, k
        real(kind=8) :: t1,t2

        ! Read adjoint specific information
        
        levtol = 0.0d0

        inquire(file=adjointfile, exist=fileExists)
        if (fileExists) then

            iunit = 16
            call opendatafile(iunit,adjointfile)
            
            read(iunit,*) adjoint_flagging
            if (.not. adjoint_flagging) return
            
            read(iunit,*) adjoint_output

            ! time period of interest:
            read(iunit,*) t1
            read(iunit,*) t2

            call set_time_window(t1, t2)

            read(iunit,*) innerprod_index
            
            read(iunit,*) totnum_adjoints
            allocate(adj_files(totnum_adjoints))

            do 20 k = 1, totnum_adjoints
                read(iunit,*) adj_files(totnum_adjoints + 1 - k)
            20  continue
            close(iunit)

            if (adjoint_flagging .and. (size(adj_files) <= 0)) then
                print *, 'Error: no adjoint output files found.'
                stop
            endif

            ! Allocate space for the number of needed binary output files
            allocate(adjoints(totnum_adjoints))

            do 50 k = 1, totnum_adjoints
                ! Load binary output files
                call reload(adj_files(k),k)

            50 continue

        else
            adjoint_flagging = .false.
        endif

    end subroutine read_adjoint_data

! ========================================================================
!  set_time_window(t1,t2)
! ========================================================================
    subroutine set_time_window(t1, t2)

        implicit none

        ! Function Arguments
        real(kind=8), intent(in) :: t1, t2

        trange_final = t2
        trange_start = t1

    end subroutine set_time_window

! ========================================================================
!  calculate_tol(lcheck)
! ========================================================================
    subroutine calculate_tol(lcheck)

        use amr_module
        implicit none

        ! Function Arguments
        integer, intent(in) :: lcheck
        real(kind=8) :: errtotal, cutoff
        real(kind=8) :: dt,hx,hy
        integer :: celln, sorted(numcells(lcheck)/2)

        levtol(lcheck) = NEEDS_TO_BE_SET

        dt = possk(lcheck)
        hx = hxposs(lcheck)
        hy = hyposs(lcheck)

        ! Setting our goal for the maximum amount of error
        ! for this level
        cutoff = tol*(dt/(hx*hy))/(tfinal-t0)
        cutoff = cutoff/mxnest

        ! Sorting errors
        call qsortr(sorted, numcells(lcheck)/2, errors)

        errtotal = 0
        do celln = 1, numcells(lcheck)/2
            errtotal = errtotal + errors(sorted(celln))
            if (errtotal .ge. cutoff) then
                levtol(lcheck) = errors(sorted(celln-1))
                EXIT
            endif
        end do

    end subroutine calculate_tol

! ========================================================================
!  reload(adjfile, k)
!  Note: This assumes that the binary output format was used
! ========================================================================

    subroutine reload(adjfile, k)

        implicit double precision (a-h,o-z)

        integer, intent(in) :: k
        integer :: mptr, level, ladjfile, mptr_notused
        integer :: mitot, mjtot, i1, i2
        integer :: i,j, ivar, z, loc
        integer :: allocsize, new_size
        character(len=*) :: adjfile
        logical foundFile, initial

        real(kind=8), allocatable, target, dimension(:) :: new_storage

        iadd(ivar,i,j)  = adjoints(k)%loc(mptr) &
            + ivar - 1 + adjoints(k)%meqn*((j-1)*mitot+i-1)

       ! Checking to see if fort.t file exists
        ladjfile = len(trim(adjfile))
        adjfile(ladjfile-4:ladjfile-4) = 't'
        !write(6,*) 'Attempting to reload data '
        !write(6,*) '  fort.t* file: ',trim(adjfile)
        inquire(file=trim(adjfile),exist=foundFile)
        if (.not. foundFile) then
            write(*,*)" Did not find fort.t* file!"
            stop
        endif
        open(9,file=trim(adjfile),status='old',form='formatted')
        rewind 9

        ! Reading from fort.t file
        read(9, *) adjoints(k)%time
        read(9, *) adjoints(k)%meqn
        read(9, *) adjoints(k)%ngrids
        read(9, *) adjoints(k)%naux
        read(9, *) adjoints(k)%ndim
        read(9, *) adjoints(k)%nghost

        close(9)
        write(*,*) 'Loading adjoint data at time t = ', adjoints(k)%time

       ! Allocating memory for alloc array
        allocsize = 4000000
        allocate(adjoints(k)%alloc(allocsize))

      ! Checking to see if fort.q file exists
        adjfile(ladjfile-4:ladjfile-4) = 'q'
        !write(6,*) 'Attempting to reload data '
        !write(6,*) '  fort.q* file: ',trim(adjfile)
        inquire(file=trim(adjfile),exist=foundFile)
        if (.not. foundFile) then
            write(*,*)" Did not find fort.q* file!"
            stop
        endif
        open(10,file=trim(adjfile),status='old',form='formatted')
        rewind 10

       ! Checking to see if fort.b file exists
        adjfile(ladjfile-4:ladjfile-4) = 'b'
        !write(6,*) 'Attempting to reload data '
        !write(6,*) '  fort.b* file: ',trim(adjfile)
        inquire(file=trim(adjfile),exist=foundFile)
        if (.not. foundFile) then
            write(*,*)" Did not find fort.b* file!"
            stop
        endif
        write(*,*) '   from file ', trim(adjfile)
        write(*,*) ' '
        open(20,file=trim(adjfile),status='unknown',access='stream')
        rewind 20

       ! Allocating size for grid information arrays
       allocate(adjoints(k)%xlowvals(adjoints(k)%ngrids))
       allocate(adjoints(k)%ylowvals(adjoints(k)%ngrids))
       allocate(adjoints(k)%ncellsx(adjoints(k)%ngrids))
       allocate(adjoints(k)%ncellsy(adjoints(k)%ngrids))
       allocate(adjoints(k)%loc(adjoints(k)%ngrids))
       allocate(adjoints(k)%gridlevel(adjoints(k)%ngrids))
       allocate(adjoints(k)%gridpointer(adjoints(k)%ngrids))

       ! Initializing all levels to zero
       adjoints(k)%gridlevel(:) = 0

       ! Reading from fort.q* file and fort.b* files
        loc = 1
        do z = 1, adjoints(k)%ngrids
            read(10,*) mptr_notused
            adjoints(k)%gridpointer(z) = z
            mptr = z

            read(10,*) level
            adjoints(k)%gridlevel(mptr) = level

            read(10,*) adjoints(k)%ncellsx(mptr)
            read(10,*) adjoints(k)%ncellsy(mptr)
            read(10,*) adjoints(k)%xlowvals(mptr)
            read(10,*) adjoints(k)%ylowvals(mptr)
            read(10,*) adjoints(k)%hxposs(level)
            read(10,*) adjoints(k)%hyposs(level)
            read(10,*)

            mitot = adjoints(k)%ncellsx(mptr) + 2*adjoints(k)%nghost
            mjtot = adjoints(k)%ncellsy(mptr) + 2*adjoints(k)%nghost

            adjoints(k)%loc(mptr) = loc
            loc = loc + mitot*mjtot*adjoints(k)%meqn

            ! Checking to see if the alloc array is large enough
            ! to hold the new grid
            ! If not, making the alloc array larger
            if (allocsize .lt. loc) then
                new_size = 2*allocsize
                allocate(new_storage(new_size))

                new_storage(1:allocsize) = adjoints(k)%alloc
                call move_alloc(new_storage,adjoints(k)%alloc)
                allocsize = new_size
            endif

            ! This is the bulk of the reading
            i1 = iadd(1,1,1)
            i2 = iadd(adjoints(k)%meqn,mitot,mjtot)
            read(20) adjoints(k)%alloc(i1:i2)

        enddo

        close(10)
        close(20)

        adjoints(k)%lfine = maxval(adjoints(k)%gridlevel)

    end subroutine reload

! ========================================================================
!  select_snapshots(time,mask_selecta)
! ========================================================================
    subroutine select_snapshots(time,mask_selecta)

       implicit none

       logical, intent(inout) :: mask_selecta(totnum_adjoints)
       real(kind=8), intent(in) :: time

!      Local variables
       logical :: adjoints_found
       integer :: r

!       Pick adjoint snapshots to consider when flagging
        mask_selecta = .false.
        adjoints_found = .false.

        do r=1,totnum_adjoints
            if ((time+adjoints(r)%time) >= trange_start .and. &
              (time+adjoints(r)%time) <= trange_final) then
                mask_selecta(r) = .true.
                adjoints_found = .true.
            endif
        enddo

        if(.not. adjoints_found) then
            write(*,*) "Error: no adjoint snapshots ", &
                "found in time range."
            write(*,*) "Consider increasing time rage of interest, ", &
                "or adding more snapshots."
        endif

        do r=1,totnum_adjoints-1
            if((.not. mask_selecta(r)) .and. &
              (mask_selecta(r+1))) then
                mask_selecta(r) = .true.
                exit
            endif
        enddo

        do r=totnum_adjoints,2,-1
            if((.not. mask_selecta(r)) .and. &
              (mask_selecta(r-1))) then
                mask_selecta(r) = .true.
                exit
            endif
        enddo

    end subroutine select_snapshots

end module adjoint_module