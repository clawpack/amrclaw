!
!> Check each grid, starting with **mptr1** (either newstl or lstart)
!! to see that it has no more than max1d points in either dimensions.
!! needed so that scratch array space in stepgrid not exceeded.
!!
!! Also check for too small grids - but has never happened.
!
! --------------------------------------------------
!
subroutine birect(mptr1)
    !
    use amr_module
    implicit real(CLAW_REAL) (a-h,o-z)
    integer :: total_cells


    !
    ! :::::::::::::  BIRECT :::::::::::::::::::::::::::::::::::::::
    ! check each grid, starting with mptr1 (either newstl or lstart)
    ! to see that it has no more than max1d points in either dimensions.
    ! needed so that scratch array space in stepgrid not exceeded.
    !
    ! also check for too small grids - but has never happened.
    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    mptr  = mptr1
    level = node(nestlevel,mptr)
    hx    = hxposs(level)
    hy    = hyposs(level)
    !
    10    continue
    cxlo    = rnode(cornxlo,mptr)
    cxhi    = rnode(cornxhi,mptr)
    cylo    = rnode(cornylo,mptr)
    cyhi    = rnode(cornyhi,mptr)
    nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
    ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
    minsize = 2*nghost
    !
    ! check number of rows first - if too many, bisect grid with vertical
    ! line down the middle. make sure new grid corners are anchored
    ! on coarse grid point. make sure if bisecting coarse grid that
    ! new grids have even number of points
    !
    total_cells = (nx+2*nghost)*(ny+2*nghost)
    if ( &
         (nx + 2*nghost .gt. max1d ) &
         .and. (total_cells > max1d*max1d/2) &
       ) then

        nxl    = nx/2
        if (level .gt. 1) then 
            lratio = intratx(level-1)
        else
            lratio = 2
        endif
        nxl = (nxl/lratio)*lratio 
        nxr    = nx - nxl 
        cxmid  = cxlo + nxl*hx

        mptrnx = nodget()
        node(levelptr,mptrnx) = node(levelptr,mptr)
        node(levelptr,mptr)   = mptrnx

        rnode(cornxhi,mptr) = cxmid
        node(ndihi,mptrnx)  = node(ndihi,mptr)
        node(ndihi,mptr)    = node(ndilo,mptr) + nxl - 1
        node(ndilo,mptrnx)  = node(ndihi,mptr) + 1
        node(ndjhi,mptrnx)  = node(ndjhi,mptr)
        node(ndjlo,mptrnx)  = node(ndjlo,mptr)

        rnode(cornxlo,mptrnx)    = cxmid
        rnode(cornylo,mptrnx)    = cylo
        rnode(cornyhi,mptrnx)    = cyhi
        rnode(cornxhi,mptrnx)    = cxhi
        rnode(timemult,mptrnx)   = rnode(timemult,mptr)
        node(nestlevel,mptrnx)   = node(nestlevel,mptr)

        go to 10
        !
        ! check number of columns next - if too many, bisect grid with horizontal
        ! line down the middle
        !
    else if ( &
              (ny + 2*nghost .gt. max1d) &
              .and. (total_cells > max1d*max1d/2) &
            ) then

        nyl    = ny/2
        if (level .gt. 1) then 
            lratio = intraty(level-1)
        else
            lratio = 2
        endif
        nyl = (nyl/lratio)*lratio
        nyr    = ny - nyl 
        cymid  =  cylo + nyl*hy

        mptrnx = nodget()
        node(levelptr,mptrnx) = node(levelptr,mptr)
        node(levelptr,mptr)   = mptrnx

        rnode(cornyhi,mptr)   = cymid

        node(ndjhi,mptrnx) = node(ndjhi,mptr)
        node(ndjhi,mptr)   = node(ndjlo,mptr) + nyl - 1
        node(ndjlo,mptrnx) = node(ndjhi,mptr) + 1
        node(ndihi,mptrnx) = node(ndihi,mptr)
        node(ndilo,mptrnx) = node(ndilo,mptr)

        rnode(cornxlo,mptrnx)   = cxlo
        rnode(cornylo,mptrnx)   = cymid
        rnode(cornyhi,mptrnx)   = cyhi
        rnode(cornxhi,mptrnx)   = cxhi
        node(nestlevel,mptrnx)  = node(nestlevel,mptr)
        rnode(timemult,mptrnx)  = rnode(timemult,mptr)
        go to 10
        !
        !  grid ok - check the next
        !
    else
        mptr = node(levelptr,mptr)
        if (mptr.ne.0) go to 10

    endif

    return
end
