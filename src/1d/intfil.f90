!
! ::::::::::::::::::::::: INTFIL ::::::::::::::::::::::::::::::::;
!  INTFIL: interpolates values for a patch at the specified level and
!  location, using values from grids at LEVEL and coarser, if nec.
!
!  take the intersection of a grid patch with corners at ilo,ihi
!  and all grids mptr at LEVEL.  If there is a non-null intersection
!  copy the solution vaues from mptr (at TIME) into VAL array.
!  assumes patch at same levels do straight copy, not skipping
!  every intrat or doing any interpolation here,
!  assume called in correct order of levels, so that when copying
!  is ok to overwrite.
!
!  N.B.: there are no dummy points around patch, since
!        this is not an official "grid" but a "patch".
!
!  used array marks when point filled. filpatch checks if points left over
!  after intersections at specified level.
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;
!
subroutine intfil(val,mi,time,flaguse,nrowst,ilo,ihi,level,nvar,naux,msrc)

    use amr_module, only: possk, mxnest, iregsz, nghost, outunit, alloc
    use amr_module, only: node, lstart, store1, store2, levelptr, timemult,gridNbor
    use amr_module, only: rnode, ndilo, ndihi, nextfree
    use amr_module, only: bndListNum, bndListSt
    use amr_module, only: listStart, listOfGrids, bndList, numgrids

    implicit none

    ! Input
    integer, intent(in) :: mi, nrowst, ilo, ihi, level, nvar, naux,msrc
    real(kind=8), intent(in) :: time

    ! In/Out
    integer(kind=1), intent(in out) :: flaguse(ilo:ihi)
    real(kind=8), intent(in out) :: val(nvar,mi)

    ! Locals
    integer :: imlo, imhi, nx, mitot
    integer :: ixlo, ixhi, locold, locnew, nextSpot
    integer :: icount, bndNum, bndLoc, levSt
    integer :: ivar, i, mptr, mstart, loc, numg
    real(kind=8) :: dt, alphac, alphai
    logical :: t_interpolate

    integer :: patch_line(2)

    real(kind=8), parameter :: t_epsilon = 1.0d-4

    ! Formats for error statements
    character(len=*), parameter :: missing_error = &
            "(' time wanted ',e15.7,' not available from grid ',i4,'level',i4)"
    character(len=*), parameter :: time_error = &
            "(' trying to interpolate from previous time values ',/," // &
            "' for a patch with corners ilo,ihi,jlo,jhi:'" // &
            ",/,2x,4i10,/," // &
            "' from source grid ',i4,' at level ',i4,/," // &
            "' time wanted ',e24.16,' source time is ',e24.16,/," // &
            "' alphai, t_epsilon ',2e24.16)"

    patch_line = [ilo,ihi]

    ! Note that we need a non-dimensionalized t epspatch_line(1)n as it was a problem
    ! in tsunami tests ( instead of t_epsilon   = dt / 10.d0 )
    
    ! Time step at this level
    dt = possk(level)
      
    ! Initialize the flagging where we set things
    flaguse = 0

    ! Loop through all grids at this level, initialize to first
!    mptr = lstart(level)
!    do while (mptr /= 0)
     if (msrc .eq. -1) then
         numg = numgrids(level)
         levSt = listStart(level)
     else
         bndLoc = node(bndListSt,msrc)  ! index of first grid in bndList
         bndNum = node(bndListNum,msrc)
         nextSpot = node(bndListSt, msrc) ! initialize
         numg = bndNum
     endif

     do icount = 1, numg

         if (msrc .eq. -1) then
            mptr = listOfGrids(levSt+icount-1)
         else
            mptr = bndList(nextSpot,gridNbor)
         endif

        ! Check if grid mptr and patch intersect
        imlo = node(ndilo, mptr)
        imhi = node(ndihi, mptr)

        nx = node(ndihi,mptr) - node(ndilo,mptr) + 1
        
        mitot = nx + 2 * nghost

        ixlo = max(imlo,patch_line(1))
        ixhi = min(imhi,patch_line(2))

        ! Check to see if grid and patch interesect, if not continue to next
        ! grid in the list
        if (ixlo <= ixhi) then

            ! grids intersect. figure out what time to use.
            ! alphai = 1 for new time; 0 for old time
            alphac = (rnode(timemult,mptr) - time)/dt
            alphai = 1.d0 - alphac

            if ((alphai < -t_epsilon) .or. (alphai > 1.d0 + t_epsilon)) then
                write(outunit,missing_error) time, mptr, level
                print missing_error, time, mptr, level
                write(outunit,'(A,E24.16)') 'Line 80', dt
                write(outunit,time_error) patch_line,mptr,level,time,rnode(timemult,mptr),alphai,t_epsilon
                print time_error, patch_line,mptr,level,time,rnode(timemult,mptr),alphai,t_epsilon
                call outtre(mstart,.false.,nvar,naux)
                stop
            endif

            ! Check if we should interpolate in time
            t_interpolate = .false.
            if (abs(alphai - 1.d0) < t_epsilon) then
                loc = node(store1,mptr)
            else if (dabs(alphai) .lt. t_epsilon) then
                loc = node(store2,mptr)
                if (level == mxnest) then
                    write(outunit,'(A,E24.16)') 'Line 95', dt
                    write(outunit,time_error) patch_line,mptr,level,time, &
                                              rnode(timemult,mptr),alphai,t_epsilon
                    write(*,time_error) patch_line,mptr,level,time, &
                                        rnode(timemult,mptr),alphai,t_epsilon
                    stop
                endif
            else
                locold  = node(store2,mptr)
                locnew  = node(store1,mptr)
                t_interpolate = .true.

                ! If we are at the maximum level nesting, abort
                if (level == mxnest) then
                    write(outunit,'(A,E24.16)') 'Line 107',dt
                    write(outunit,time_error) patch_line,mptr,level,time,rnode(timemult,mptr),alphai,t_epsilon
                    print time_error, patch_line,mptr,level,time,rnode(timemult,mptr),alphai,t_epsilon
                    stop
                endif
            endif

            ! Actual interpolation
            if (.not. t_interpolate) then
                ! No time interp. copy the solution values
                do ivar = 1, nvar
                    do i = ixlo, ixhi
                        val(ivar,i-patch_line(1)+nrowst) = &
                            alloc(iadd(ivar,i-imlo+nghost+1))
                        flaguse(i) = 1
                    end do
                end do
            else
                ! Linear interpolation in time
                do ivar = 1, nvar
                    do i = ixlo, ixhi
                        val(ivar,i-patch_line(1)+nrowst) = &
                            alloc(iadnew(ivar,i-imlo+nghost+1))*alphai + &
                            alloc(iadold(ivar,i-imlo+nghost+1))*alphac
                        flaguse(i) = 1
                    end do
                end do
            endif

        endif

        ! Get next grid
!        mptr = node(levelptr, mptr)
         if (msrc .ne. -1) then
            nextSpot = bndList(nextSpot,nextfree)
         endif

    end do

    ! Set used array points which intersect domain boundary to be equal to 1, 
    ! they will be set elsewhere, namely boundary conditions and periodic
    ! domains

    if (patch_line(1) < 0) then
        flaguse(patch_line(1):min(-1,nrowst + ihi - patch_line(1))) = 1
    endif

    if (ihi >= iregsz(level)) then
        flaguse(max(iregsz(level),patch_line(1)):ihi) = 1
    endif

contains

    integer pure function iadd(ivar,i)
        implicit none
        integer, intent(in) :: ivar, i
        iadd = loc + ivar-1 + nvar*(i-1)
    end function iadd

    integer pure function iadnew(ivar,i)
        implicit none
        integer, intent(in) :: ivar, i
        iadnew = locnew + ivar-1 + nvar*(i-1)
    end function iadnew

    integer pure function iadold(ivar,i)
        implicit none
        integer, intent(in) :: ivar, i
        iadold = locold + ivar-1 + nvar*(i-1)
    end function iadold

end subroutine intfil
