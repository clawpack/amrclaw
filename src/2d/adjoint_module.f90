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

            adjoint_flagging = .true.
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

            ! Allocate space for the number of needed binary output files
            allocate(adjoints(totnum_adjoints))

            do 50 k = 1, totnum_adjoints
                ! Load binary output files
                call reload(adj_files(k),k)

            50 continue

        else
            print *, 'Adjoint.data file does not exist.'
            print *, 'If you are using the adjoint method, this is an error.'
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
        write(6,*) 'Attempting to reload data '
        write(6,*) '  fort.t* file: ',trim(adjfile)
        inquire(file=trim(adjfile),exist=foundFile)
        if (.not. foundFile) then
            write(*,*)" Did not find fort.t* file!"
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

       ! Allocating memory for alloc array
        allocsize = 4000000
        allocate(adjoints(k)%alloc(allocsize))

      ! Checking to see if fort.q file exists
        adjfile(ladjfile-4:ladjfile-4) = 'q'
        write(6,*) 'Attempting to reload data '
        write(6,*) '  fort.q* file: ',trim(adjfile)
        inquire(file=trim(adjfile),exist=foundFile)
        if (.not. foundFile) then
            write(*,*)" Did not find fort.q* file!"
            stop
        endif
        open(10,file=trim(adjfile),status='old',form='formatted')
        rewind 10

       ! Checking to see if fort.b file exists
        adjfile(ladjfile-4:ladjfile-4) = 'b'
        write(6,*) 'Attempting to reload data '
        write(6,*) '  fort.b* file: ',trim(adjfile)
        inquire(file=trim(adjfile),exist=foundFile)
        if (.not. foundFile) then
            write(*,*)" Did not find fort.b* file!"
            stop
        endif
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
            read(10,"(i6)") mptr_notused
            adjoints(k)%gridpointer(z) = z
            mptr = z

            read(10,"(i6)") level
            adjoints(k)%gridlevel(mptr) = level

            read(10,"(i6)") adjoints(k)%ncellsx(mptr)
            read(10,"(i6)") adjoints(k)%ncellsy(mptr)
            read(10,"(e26.16)") adjoints(k)%xlowvals(mptr)
            read(10,"(e26.16)") adjoints(k)%ylowvals(mptr)
            read(10,"(e26.16)") adjoints(k)%hxposs(level)
            read(10,"(e26.16)") adjoints(k)%hyposs(level)
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

! ========================================================================
!  Routine to calculate inner product
! ========================================================================

    subroutine calculate_innerproduct(q,k,mx_f,my_f,xlower_f, &
               ylower_f,dx_f,dy_f,meqn_f,mbc_f, innerprod)

        implicit none

        real(kind=8), intent(in) :: xlower_f,ylower_f,dx_f,dy_f
        integer :: k,mx_f,my_f,meqn_f,mbc_f
        real(kind=8), intent(in) :: q(meqn_f,1-mbc_f:mx_f+mbc_f,1-mbc_f:my_f+mbc_f)

        integer :: mx_a, my_a, mptr_a, mbc_a
        integer :: i, j, i1, i2, j1, j2, level, loc, z
        real(kind=8) :: dx_a, xlower_a, xupper_a, xupper_f
        real(kind=8) :: dy_a, ylower_a, yupper_a, yupper_f
        real(kind=8) :: x1, x2, y1, y2

        real(kind=8), intent(inout) :: innerprod(mx_f,my_f)
        real(kind=8) :: q_innerprod(mx_f,my_f)
        logical :: mask_forward(mx_f,my_f)
        real(kind=8) :: q_interp(meqn_f,mx_f,my_f)

        logical, allocatable :: mask_adjoint(:,:)

        xupper_f = xlower_f + mx_f*dx_f
        yupper_f = ylower_f + my_f*dy_f

        ! Loop over patches in adjoint solution
        do z = 1, adjoints(k)%ngrids
            mptr_a = adjoints(k)%gridpointer(z)
            level = adjoints(k)%gridlevel(mptr_a)

            ! Number of points in x and y
            mx_a = adjoints(k)%ncellsx(mptr_a)
            my_a = adjoints(k)%ncellsy(mptr_a)

            ! Finding extreem values for grid
            xlower_a = adjoints(k)%xlowvals(mptr_a)
            dx_a = adjoints(k)%hxposs(level)
            xupper_a = xlower_a + mx_a*dx_a
            ylower_a = adjoints(k)%ylowvals(mptr_a)
            dy_a = adjoints(k)%hyposs(level)
            yupper_a = ylower_a + my_a*dy_a

            loc = adjoints(k)%loc(mptr_a)

            ! Check if adjoint patch overlaps with forward patch
            x1 = max(xlower_f,xlower_a)
            x2 = min(xupper_f,xupper_a)
            y1 = max(ylower_f,ylower_a)
            y2 = min(yupper_f,yupper_a)

            if ((x1 > x2) .or. (y1 > y2)) then
                ! Skipping interpolation if grids don't overlap
                mask_forward = .false.
                continue
            else
                mbc_a = adjoints(k)%nghost
                allocate(mask_adjoint(1-mbc_a:mx_a+mbc_a, 1-mbc_a:my_a+mbc_a))

                ! Create a mask that is .true. only in part of patch intersecting forward patch:
                i1 = max(int((x1 - xlower_a + 0.5d0*dx_a) / dx_a), 0)
                i2 = min(int((x2 - xlower_a + 0.5d0*dx_a) / dx_a) + 1, mx_a+1)
                j1 = max(int((y1 - ylower_a + 0.5d0*dy_a) / dy_a), 0)
                j2 = min(int((y2 - ylower_a + 0.5d0*dy_a) / dy_a) + 1, my_a+1)

                forall (i=1-mbc_a:mx_a+mbc_a, j=1-mbc_a:my_a+mbc_a)
                    mask_adjoint(i,j) = ((i >= i1) .and. (i <= i2) .and. &
                                       (j >= j1) .and. (j <= j2))
                end forall

                ! Create a mask that is .true. only in part of forward patch intersecting patch:

                i1 = max(int((x1 - xlower_f + 0.5d0*dx_f) / dx_f)+1, 0)
                i2 = min(int((x2 - xlower_f + 0.5d0*dx_f) / dx_f), mx_f)
                j1 = max(int((y1 - ylower_f + 0.5d0*dy_f) / dy_f)+1, 0)
                j2 = min(int((y2 - ylower_f + 0.5d0*dy_f) / dy_f), my_f)

                forall (i=1:mx_f, j=1:my_f)
                    mask_forward(i,j) = ((i >= i1) .and. (i <= i2) .and. &
                                       (j >= j1) .and. (j <= j2))
                end forall

                ! Interpolate adjoint values to q_interp
                ! Note that values in q_interp will only be set properly where
                ! mask_adjoint == .true.
                call interp_adjoint( &
                        adjoints(k)%meqn, k, q_interp, xlower_a, ylower_a, &
                        dx_a, dy_a, mx_a, my_a, xlower_f, ylower_f, &
                        dx_f, dy_f, mx_f, my_f, &
                        mask_adjoint, mptr_a, mask_forward)

                q_innerprod = 0.d0
                ! For each overlapping point, calculate inner product
                forall(i = 1:mx_f, j = 1:my_f, mask_forward(i,j))
                    q_innerprod(i,j) = abs(dot_product(q(:,i,j),q_interp(:,i,j)))
                end forall

                do i=1,mx_f
                    do j=1,my_f
                        if (q_innerprod(i,j) > innerprod(i,j)) then
                            innerprod(i,j) = q_innerprod(i,j)
                        endif
                    enddo
                enddo

                deallocate(mask_adjoint)
            endif
        enddo

    end subroutine calculate_innerproduct

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! :::::     Routine to interpolate adjoint to given x,y point
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine interp_adjoint(nvar, r, q_interp, xlower_a, ylower_a, dx_a, dy_a, &
               mx_a, my_a, xlower_f, ylower_f, dx_f, dy_f, mx_f, my_f, &
               mask_adjoint, mptr_a, mask_forward)

        implicit none

        ! Function arguments
        integer, intent(in) :: r, nvar, mx_a, my_a
        logical, intent(in) :: mask_adjoint(1-adjoints(r)%nghost:mx_a+adjoints(r)%nghost, &
        1-adjoints(r)%nghost:my_a+adjoints(r)%nghost)
        real(kind=8), intent(in) :: xlower_a, xlower_f, ylower_a, ylower_f
        integer, intent(in) :: mx_f, my_f, mptr_a
        real(kind=8), intent(in) :: dx_f, dx_a, dy_f, dy_a

        integer :: z,k, iz, jk, mitot
        integer :: ivar, i, j, iadd, iaddaux, loc
        real(kind=8) :: q_interp(nvar,mx_f,my_f), denom
        real(kind=8) :: x, xhigh_a, y, yhigh_a
        real(kind=8) :: dxz,dyk, a, b, c
        logical :: mask_forward(mx_f,my_f)

        iadd(ivar,i,j)  = loc + ivar - 1 + adjoints(r)%meqn*((j-1)*mitot+i-1)

        q_interp = 0.0
        xhigh_a  = xlower_a + mx_a*dx_a
        yhigh_a = ylower_a + my_a*dx_a
        loc    = adjoints(r)%loc(mptr_a)
        mitot = adjoints(r)%ncellsx(mptr_a) + 2*adjoints(r)%nghost

        do z = 1,mx_f
            do k = 1,my_f
                if (mask_forward(z,k)) then
                    x = xlower_f + (z - 0.5d0)*dx_f
                    y = ylower_f + (k - 0.5d0)*dy_f

                    !TODO: Why does iz and jk have an added 1 at the end?

                    iz = int((x - xlower_a + 0.5d0*dx_a) / dx_a) + 1
                    dxz = x - (xlower_a + (iz-0.5d0)*dx_a)
                    jk = int((y - ylower_a + 0.5d0*dy_a) / dy_a) + 1
                    dyk = y - (ylower_a + (jk-0.5d0)*dy_a)

                    ! Interpolate only if this cell is overlapping with grid
                    if (mask_adjoint(iz,jk)) then
                        do ivar=1,nvar

                        a = (adjoints(r)%alloc(iadd(ivar,iz+1,jk)) &
                              - adjoints(r)%alloc(iadd(ivar,iz,jk))) / dx_a
                        b = (adjoints(r)%alloc(iadd(ivar,iz,jk+1)) &
                              - adjoints(r)%alloc(iadd(ivar,iz,jk))) / dy_a
                        c = (adjoints(r)%alloc(iadd(ivar,iz+1,jk+1)) &
                              + adjoints(r)%alloc(iadd(ivar,iz,jk)) &
                              - adjoints(r)%alloc(iadd(ivar,iz+1,jk)) &
                              - adjoints(r)%alloc(iadd(ivar,iz,jk+1))) / (dx_a*dy_a)

                        q_interp(ivar,z,k) = &
                              adjoints(r)%alloc(iadd(ivar,iz,jk)) &
                              + a*dxz + b*dyk + c*dxz*dyk

                        enddo
                    endif
                endif
            enddo
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

    subroutine errf1a(rctfine,nvar,rctcrse,mptr,mi2tot,mj2tot, &
                      mitot,mjtot,rctflg,mibuff,mjbuff,auxfine, &
                      naux)
      use amr_module
      implicit double precision (a-h,o-z)

 
      dimension  rctfine(nvar,mitot,mjtot)
      dimension  rctcrse(nvar,mi2tot,mj2tot)
      dimension  rctflg(mibuff,mjbuff)

      dimension  aux_crse(mi2tot,mj2tot)
      dimension  err_crse(nvar,mi2tot,mj2tot)
      logical mask_selecta(totnum_adjoints)
      dimension  auxfine(naux,mitot,mjtot)
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
      ybot  = rnode(cornylo,mptr)
      hy    = hyposs(levm)
      dt    = possk(levm)
      numsp = 0
      jg = grid_num(mptr)
      nx = mitot - 2*nghost
      ny = mjtot - 2*nghost

      call select_snapshots(time,mask_selecta)

      errmax = 0.0d0
      err2   = 0.0d0
      auxfine(innerprod_index,:,:) = 0.0d0
      aux_crse = 0.0d0

      order  = dble(2**(iorder+1) - 2)
!
!     Calculating correct tol for this level
!     nxnest is the maximum number of refinement levels, from amr_module
!     --------------------
      tol_exact = tol*dt/(tfinal-t0)
      tol_exact = tol_exact/mxnest
      tol_exact = tol_exact/(numcells(levm)*hx*hy)

      if (t0+possk(levm) .eq. time) levtol(levm) = tol_exact

!
! zero out the exterior locations so they don't affect err.est.
!
      jfine = nghost+1
      do j = nghost+1, mj2tot-nghost
        yofj  = ybot + (dble(jfine) - .5d0)*hy
        ifine = nghost+1

        do i  = nghost+1, mi2tot-nghost
          xofi  = xleft + (dble(ifine) - .5d0)*hx

          do k = 1, nvar
            term1 = rctfine(k,ifine,jfine)
            term2 = rctfine(k,ifine+1,jfine)
            term3 = rctfine(k,ifine+1,jfine+1)
            term4 = rctfine(k,ifine,jfine+1)
!           # divide by (aval*order) for relative error
            aval  = (term1+term2+term3+term4)/4.d0
            est   =  dabs((aval-rctcrse(k,i,j))/ order)
            if (est .gt. errmax) errmax = est
              err2 = err2 + est*est

              err_crse(k,i,j) = est
!             retaining directionality of the wave
              err_crse(k,i,j) = sign(est,rctcrse(k,i,j))
          enddo

          ifine = ifine + 2
        enddo

        jfine = jfine + 2
      enddo

      do k = 1,totnum_adjoints
!         Consider only snapshots that are within the desired time range
          if (mask_selecta(k)) then
!             set innerproduct
              call calculate_innerproduct(err_crse,k,nx/2, &
                    ny/2,xleft,ybot,hx*2,hy*2,nvar,nghost, &
                    aux_crse(nghost+1:mi2tot-nghost,nghost+1:mj2tot-nghost))
          endif
      enddo

      do i  = nghost+1, mi2tot-nghost
        do j = nghost+1, mj2tot-nghost
          i_val = i-nghost
          j_val = j-nghost
          errors(eptr(jg)+(i_val-1)*ny/2+j_val) = aux_crse(i,j)

          rctcrse(1,i,j)  = goodpt
          if (aux_crse(i,j) .ge. levtol(levm)) then
!                    ## never set rctflg to good, since flag2refine may
!                    ## have previously set it to bad
!                    ## can only add bad pts in this routine
              rctcrse(1,i,j)  = badpt
          endif

        enddo
      enddo

      jfine   = nghost+1
      do j = nghost+1, mj2tot-nghost
        ifine   = nghost+1

        do i = nghost+1, mi2tot-nghost
          auxfine(innerprod_index,ifine,jfine) = aux_crse(i,j)
          auxfine(innerprod_index,ifine+1,jfine) = aux_crse(i,j)
          auxfine(innerprod_index,ifine,jfine+1) = aux_crse(i,j)
          auxfine(innerprod_index,ifine+1,jfine+1) = aux_crse(i,j)

          if (rctcrse(1,i,j) .ne. goodpt) then
!           ## never set rctflg to good, since flag2refine may
!           ## have previously set it to bad
!           ## can only add bad pts in this routine
            rctflg(ifine,jfine)    = badpt
            rctflg(ifine+1,jfine)  = badpt
            rctflg(ifine,jfine+1)  = badpt
            rctflg(ifine+1,jfine+1)= badpt
          endif
          ifine   = ifine + 2
        enddo

        jfine   = jfine + 2
      enddo

    end subroutine errf1a


end module adjoint_module
