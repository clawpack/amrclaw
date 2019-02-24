! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::::     Module to define and work with adjoint type
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module adjoint_module
    use amr_module

    type adjointData_type

        real(kind=8), allocatable, dimension(:) :: alloc

        real(kind=8) hxposs(maxlv)

        integer meqn, ngrids, naux, ndim, nghost, lfine
        real(kind=8) time
        real(kind=8),allocatable,dimension(:) :: xlowvals
        integer, allocatable, dimension(:) :: ncellsx,loc, &
                gridlevel,gridpointer

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

        levtol = 0.0d0

        inquire(file=adjointfile, exist=fileExists)
        if (fileExists) then

            !adjoint_flagging = .true.
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
       real(kind=8) :: dt, hx
       integer :: celln, sorted(numcells(lcheck)/2)

       dt = possk(lcheck)
       hx = hxposs(lcheck)

       ! Setting our goal for the maximum amount of error
       ! for this level
       cutoff = tol*(dt/hx)/(tfinal-t0)
       cutoff = cutoff/mxnest

       ! Sorting errors
       call qsortr(sorted, numcells(lcheck)/2, errors)

       errtotal = 0
       do celln = 1, numcells(lcheck)/2
           errtotal = errtotal + errors(sorted(celln))
           if (errtotal .ge. cutoff) then
               levtol(lcheck) = errors(sorted(celln))
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
      integer :: i, ivar, z, loc
      integer :: allocsize, new_size
      character(len=*):: adjfile
      logical foundFile, initial

      real(kind=8), allocatable, target, dimension(:) :: new_storage
      iadd(ivar,i)  = adjoints(k)%loc(mptr) &
             + ivar - 1 + adjoints(k)%meqn*(i-1)

      ! Checking to see if fort.t file exists
      ladjfile = len(trim(adjfile))
      adjfile(ladjfile-4:ladjfile-4) = 't'
      !write(6,*) 'Attempting to load adjoint data '
      !write(6,*) '  fort.t* file: ',trim(adjfile)
      inquire(file=trim(adjfile),exist=foundFile)
      if (.not. foundFile) then
          write(*,*)" Did not find fort.t* file!"
          write(*,*)"   ", trim(adjfile)
          stop
      endif
      open(9,file=trim(adjfile),status='old',form='formatted')
      rewind 9

      ! Reading from fort.t file
      read(9, "(e18.8)") adjoints(k)%time
      read(9, "(i6)") adjoints(k)%meqn
      read(9, "(i6)") adjoints(k)%ngrids
      read(9, "(i6)") adjoints(k)%naux
      read(9, "(i6)") adjoints(k)%ndim
      read(9, "(i6)") adjoints(k)%nghost

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
       allocate(adjoints(k)%ncellsx(adjoints(k)%ngrids))
       allocate(adjoints(k)%loc(adjoints(k)%ngrids))
       allocate(adjoints(k)%gridlevel(adjoints(k)%ngrids))
       allocate(adjoints(k)%gridpointer(adjoints(k)%ngrids))

       !Initializing all levels to zero
       adjoints(k)%gridlevel(:) = 0

      ! Reading from fort.q* file and fort.b* files
      loc = 1
      do z = 1, adjoints(k)%ngrids
          read(10,"(i6)") mptr_notused
          adjoints(k)%gridpointer(z) = z
          mptr = z

          read(10,"(i6)") level
          adjoints(k)%gridlevel(mptr) = level

          read(10,"(i6)") adjoints(k)%ncellsx(mptr)
          read(10,"(e26.16)") adjoints(k)%xlowvals(mptr)
          read(10,"(e26.16)") adjoints(k)%hxposs(level)
          read(10,*)

          mitot = adjoints(k)%ncellsx(mptr) + 2*adjoints(k)%nghost

          adjoints(k)%loc(mptr) = loc
          loc = loc + mitot*adjoints(k)%meqn

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
          i1 = iadd(1,1)
          i2 = iadd(adjoints(k)%meqn,mitot)
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

! ========================================================================
!  Routine to calculate inner product
! ========================================================================

    subroutine calculate_innerproduct(q,r,mx_f,xlower_f, &
                dx_f,meqn_f,mbc_f,innerprod)

        implicit none

        real(kind=8), intent(in) :: xlower_f, dx_f
        integer, intent(in) :: r, mx_f, meqn_f, mbc_f
        real(kind=8), intent(in) :: q(meqn_f,1-mbc_f:mx_f+mbc_f)

        integer :: mx_a, mptr_a
        integer :: i, i1, i2, level, loc, z
        real(kind=8) :: dx_a, xlower_a, xupper_a, xupper_f, x1, x2

        real(kind=8), intent(inout) :: innerprod(1-mbc_f:mx_f+mbc_f)
        real(kind=8) :: q_innerprod(1-mbc_f:mx_f+mbc_f)
        logical :: mask_forward(1-mbc_f:mx_f+mbc_f)
        real(kind=8) :: q_interp(adjoints(r)%meqn,1-mbc_f:mx_f+mbc_f)

        logical, allocatable :: mask_adjoint(:)


        xupper_f = xlower_f + mx_f*dx_f

        ! Loop over patches in adjoint solution
        do z = 1, adjoints(r)%ngrids
            mptr_a = adjoints(r)%gridpointer(z)
            level = adjoints(r)%gridlevel(mptr_a)

            ! Number of points in x (nx)
            mx_a = adjoints(r)%ncellsx(mptr_a)

            ! Finding x extreem values for grid
            xlower_a = adjoints(r)%xlowvals(mptr_a)
            dx_a = adjoints(r)%hxposs(level)
            xupper_a = xlower_a + mx_a*dx_a

            loc = adjoints(r)%loc(mptr_a)

            ! Check if adjoint patch overlaps with forward patch
            x1 = max(xlower_f,xlower_a)
            x2 = min(xupper_f,xupper_a)

            if (x1 > x2) then
                ! Skipping interpolation if grids don't overlap
                mask_forward = .false.
                continue
            else

                allocate(mask_adjoint(1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost))

                ! Create a mask that is .true. only in part of patch intersecting forward patch:
                i1 = max(int((x1 - xlower_a + 0.5d0*dx_a) / dx_a), 0)
                i2 = min(int((x2 - xlower_a + 0.5d0*dx_a) / dx_a) + 1, mx_a+1)

                forall (i=1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost)
                    mask_adjoint(i) = ((i >= i1) .and. (i <= i2))
                end forall

                ! Create a mask that is .true. only in part of forward patch intersecting patch:

                i1 = max(int((x1 - xlower_f + 0.5d0*dx_f) / dx_f)+1, 0)
                i2 = min(int((x2 - xlower_f + 0.5d0*dx_f) / dx_f), mx_f)

                do i=1-mbc_f,mx_f+mbc_f
                    mask_forward(i) = ((i >= i1) .and. (i <= i2))
                enddo

                ! Interpolate adjoint values to q_interp
                ! Note that values in q_interp will only be set properly where 
                ! mask_adjoint == .true.
                call interp_adjoint( &
                        adjoints(r)%meqn, r, q_interp, &
                        xlower_a, dx_a, mx_a, xlower_f, dx_f, mx_f, &
                        mask_adjoint, mptr_a, mask_forward, mbc_f)

                q_innerprod = 0.d0
                ! For each overlapping point, calculate inner product
                forall(i = 1-mbc_f:mx_f+mbc_f, mask_forward(i))
                    q_innerprod(i) = abs(dot_product(q(:,i),q_interp(:,i)))
                end forall

                do i=1-mbc_f,mx_f+mbc_f
                    if (q_innerprod(i) > innerprod(i)) then
                        innerprod(i) = q_innerprod(i)
                    endif
                enddo

                deallocate(mask_adjoint)
            endif
        enddo

    end subroutine calculate_innerproduct

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given an x point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine interp_adjoint(nvar, r, q_interp, xlower_a, dx_a, &
            mx_a, xlower_f, dx_f, mx_f, mask_adjoint, mptr_a, &
            mask_forward, mbc_f)

        implicit none

        ! Function arguments
        integer, intent(in) :: r, nvar, mx_a
        logical, intent(in) :: mask_adjoint(1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost)
        real(kind=8), intent(in) :: xlower_a, xlower_f
        integer, intent(in) :: mx_f, mptr_a, mbc_f
        real(kind=8), intent(in) :: dx_f, dx_a
        logical, intent(in) :: mask_forward(1-mbc_f:mx_f+mbc_f)

        integer :: z,level, iz(1-mbc_f:mx_f+mbc_f), &
                      ivar,i, iadd, iaddaux, loc
        real(kind=8) :: q_interp(nvar,1-mbc_f:mx_f+mbc_f), denom
        real(kind=8) :: x, xhigh_a,x_f
        real(kind=8) :: dxz(1-mbc_f:mx_f+mbc_f), a

        iadd(ivar,i)  = loc + ivar - 1 + adjoints(r)%meqn*(i-1)

        q_interp = 0.0
        xhigh_a  = xlower_a + mx_a*dx_a
        loc    = adjoints(r)%loc(mptr_a)

        do z=1-mbc_f,mx_f+mbc_f
            x = xlower_f + (z - 0.5d0)*dx_f

            iz(z) = int((x - xlower_a + 0.5d0*dx_a) / dx_a) + 1
            dxz(z) = x - (xlower_a + (iz(z)-0.5d0)*dx_a)
        enddo

        do z = 1-mbc_f, mx_f+mbc_f
            x = xlower_f + (z - 0.5d0)*dx_f

            if (mask_forward(z)) then
            ! Interpolate only if this cell is overlapping with grid
                if (mask_adjoint(iz(z))) then

                    do ivar=1,nvar
                        a = (adjoints(r)%alloc(iadd(ivar,iz(z)+1)) &
                            - adjoints(r)%alloc(iadd(ivar,iz(z)))) / dx_a
                        q_interp(ivar,z) = &
                            adjoints(r)%alloc(iadd(ivar,iz(z))) + a*dxz(z)
                    enddo
                endif
            endif
        enddo

    end subroutine interp_adjoint

! ========================================================================
!  Routine to compute error in forward solution for the adjoint method
!  Editted from errf1.f
!  Differences:
!   (i) Computes error for all of the terms in q
!   (ii) Needs the aux array passed in, because the inner product is saved
!        to the aux array (note, this changes the call sequence, which is
!        why this is not included into the AMRClaw errf1.f file. If plotting
!        of the inner product is never wanted, this could be simplified
!        and included into the AMRClaw errf1.f file.
! ========================================================================
   subroutine errf1a(rctfine,nvar,rctcrse,mptr,mi2tot, &
                      mitot,rctflg,mibuff,auxfine,naux)

      implicit double precision (a-h,o-z)
 
      dimension  rctfine(nvar,mitot)
      dimension  rctcrse(nvar,mi2tot)
      dimension  rctflg(mibuff)
      dimension  aux_crse(mi2tot)
      dimension  err_crse(nvar,mi2tot)
      dimension  auxfine(naux,mitot)
      logical mask_selecta(totnum_adjoints)

!
!
! :::::::::::::::::::::: Modified from ERRF1 :::::::::::::::::::::::::
!
!  Richardson error estimator:  Used when flag_richardson is .true.
!  Compare error estimates in rctfine, rctcrse,
!  A point is flagged if the inner product between the
!  error estimate and the adjoint solution is greater than tol
!  later we check if its in a region where its allowed to be flagged
!  or alternatively required.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!

      time  = rnode(timemult, mptr)
      xleft = rnode(cornxlo,mptr)
      levm  = node(nestlevel, mptr)
      hx    = hxposs(levm)
      dt    = possk(levm)
      numsp = 0
      jg = grid_num(mptr)
      nx = mitot - 2*nghost

      call select_snapshots(time,mask_selecta)

      errmax = 0.0d0
      err2   = 0.0d0
      err_crse = 0.0d0
      auxfine(innerprod_index,:) = 0.0d0
      aux_crse = 0.0d0

      order  = dble(2**(iorder+1) - 2)
!
!     Calculating correct tol for this level
!     nxnest is the maximum number of refinement levels, from amr_module
!     --------------------
      tol_exact = tol*dt/(tfinal-t0)
      tol_exact = tol_exact/mxnest
      tol_exact = tol_exact/(numcells(levm)*hx)

      if (t0+possk(levm) .eq. time) levtol(levm) = tol_exact

!
! zero out the exterior locations so they don't affect err.est.
!
      ifine = nghost+1
      do i  = nghost+1, mi2tot-nghost

! Only check errors if flag hasn't been set yet.
! If flag == DONTFLAG then refinement is forbidden by a region,
! if flag == DOFLAG checking is not needed

! Note: here rctcrse is being used as a temporary flag
! the fine grid amrflags array is stored in rctflg, and will be
! updated based on rctcrse at the end of this routine
        if(rctflg(ifine) == UNSET &
              .or. rctflg(ifine+1) == UNSET) then
          do k = 1,nvar
              xofi  = xleft + (dble(ifine) - .5d0)*hx
              term1 = rctfine(k,ifine)
              term2 = rctfine(k,ifine+1)
!             # divide by (aval*order) for relative error
              aval  = (term1+term2)/2.d0
              est   =  dabs((aval-rctcrse(k,i))/ order)
              if (est .gt. errmax) errmax = est
              err2 = err2 + est*est

              err_crse(k,i) = est
!             retaining directionality of the wave
              err_crse(k,i) = sign(est,rctcrse(k,i))
          enddo
        else
            err_crse(:,i) = 0.d0
        endif

        ifine = ifine + 2
      enddo

      do 12 k = 1,totnum_adjoints
!         ! Consider only snapshots that are within the desired time range
          if (mask_selecta(k)) then

!             set innerproduct
              call calculate_innerproduct(err_crse,k,nx/2, &
                   xleft,hx*2,nvar,nghost,aux_crse)
          endif
 12   continue

      do i  = nghost+1, mi2tot-nghost
          errors(eptr(jg)+i-nghost) = aux_crse(i)

          rctcrse(1,i)  = DONTFLAG
          if (aux_crse(i) .ge. levtol(levm)) then
!                    ## never set rctflg to good, since flag2refine or
!                    ## flagregions may have previously set it to bad
!                    ## can only add bad pts in this routine
              rctcrse(1,i)  = DOFLAG
          endif

      enddo

      ifine   = nghost+1
      do i = nghost+1, mi2tot-nghost
         auxfine(innerprod_index,ifine) = aux_crse(i)
         auxfine(innerprod_index,ifine+1) = aux_crse(i)
         if (rctcrse(1,i) .eq. DOFLAG) then
!           ## never set rctflg to good, since flag2refine may
!           ## have previously set it to bad
!           ## can only add bad pts in this routine
            rctflg(ifine)    = DOFLAG
            rctflg(ifine+1)  = DOFLAG
        endif
        ifine   = ifine + 2
      enddo

    end subroutine errf1a


end module adjoint_module
