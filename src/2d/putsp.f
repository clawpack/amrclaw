c
c ----------------------------------------------------------
!> Reclaim list space in nodes cfluxptr and ffluxptr for all grids at level
!! **level**
c
      subroutine putsp(lbase,level,nvar,naux)
c
      use amr_module
      implicit double precision (a-h,o-z)

c
c ::::::::::::::::::::::::::::::: PUTSP :::::::::::::::::::::::::
c
c reclaim list space in nodes cfluxptr and ffluxptr for grids at level
c
c first compute max. space allocated in node cfluxptr.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      if (level .eq. lfine) go to 30
c
      ! new way - all coarse space allocated at once so
      ! reclaim it all. 
      call reclam(listspStart(level),5*listsp(level))
      ! still need to zero out coarse flux ptrs in case no
      ! new grids at next level were created (otherwise
      ! it would be re-created and assigned in prepc)
      mptr = lstart(level)
 20      node(cfluxptr,mptr) = 0
         mptr = node(levelptr,mptr)
         if (mptr .ne. 0) go to 20
      
c
 30   if (level .eq. lbase) go to 99
      mptr = lstart(level)
 40       nx    = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny    = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          ikeep = nx/intratx(level-1)
          jkeep = ny/intraty(level-1)
          lenbc = 2*(ikeep+jkeep)
c         twice perimeter since saving plus or minus fluxes 
c         plus coarse solution storage
          call reclam(node(ffluxptr,mptr), 2*nvar*lenbc+naux*lenbc)
          mptr  = node(levelptr,mptr)
          if (mptr .ne. 0) go to 40
c
 99    return
       end
