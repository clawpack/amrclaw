!
! ----------------------------------------------------------
!> Reclaim list space in nodes cfluxptr and ffluxptr for all grids at level
!! **level**
!
subroutine putsp(lbase,level,nvar,naux)
    !
    use amr_module
    implicit double precision (a-h,o-z)

    !
    ! ::::::::::::::::::::::::::::::: PUTSP :::::::::::::::::::::::::
    !
    ! reclaim list space in nodes cfluxptr and ffluxptr for grids at level
    !
    ! first compute max. space allocated in node cfluxptr.
    !
    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    if (level .ne. lfine) then
        mptr  = lstart(level)
        do while (mptr .ne. 0)
            call reclam(node(cfluxptr,mptr), 5*listsp(level))
            node(cfluxptr,mptr) = 0
            mptr  = node(levelptr,mptr)
        enddo
    endif
    !
    if (level .ne. lbase) then
        mptr = lstart(level)
        do while (mptr .ne. 0)
            nx    = node(ndihi,mptr) - node(ndilo,mptr) + 1
            ny    = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
            ikeep = nx/intratx(level-1)
            jkeep = ny/intraty(level-1)
            lenbc = 2*(ikeep+jkeep)
            !         twice perimeter since saving plus or minus fluxes 
            !         plus coarse solution storage
            call reclam(node(ffluxptr,mptr), 2*nvar*lenbc+naux*lenbc)
            mptr  = node(levelptr,mptr)
        enddo
    endif
    !
    return
end subroutine putsp
