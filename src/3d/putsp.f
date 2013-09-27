c
c ----------------------------------------------------------
c
      subroutine putsp(lbase,level,nvar,naux)
c
      use amr_module
      implicit double precision (a-h,o-z)

      parameter(numbcs=6)
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
      mptr  = lstart(level)
 20      call reclam(node(cfluxptr,mptr), numbcs*listsp(level))
         node(cfluxptr,mptr) = 0
      mptr  = node(levelptr,mptr)
      if (mptr .ne. 0) go to 20
c
 30    if (level .eq. lbase) go to 99
      mptr = lstart(level)
 40       nx    = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ny    = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
          nz    = node(ndkhi,mptr) - node(ndklo,mptr) + 1
          ikeep = nx/intratx(level-1)
          jkeep = ny/intraty(level-1)
          kkeep = nz/intratz(level-1)
          lenbc = 2*(ikeep*jkeep + jkeep*kkeep + kkeep*ikeep)
c         twice perimeter since saving plus or minus fluxes 
c         plus coarse solution storage
          call reclam(node(ffluxptr,mptr), 2*nvar*lenbc+naux*lenbc)
          mptr  = node(levelptr,mptr)
          if (mptr .ne. 0) go to 40
c
 99    return
       end
