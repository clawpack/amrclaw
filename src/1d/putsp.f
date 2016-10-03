c
c ----------------------------------------------------------
c
      subroutine putsp(lbase,level,nvar,naux)
c
      use amr_module
      implicit double precision (a-h,o-z)

      parameter(numbcs=2)

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
c BND:TODO
c This used to be "5" for 2d code, and is "numbcs" (6) for 3d code.
c So, should it be 2 or 3?
 20      call reclam(node(cfluxptr,mptr), numbcs*listsp(level))
         node(cfluxptr,mptr) = 0
      mptr  = node(levelptr,mptr)
      if (mptr .ne. 0) go to 20
c
 30    if (level .eq. lbase) go to 99
      mptr = lstart(level)
 40       nx    = node(ndihi,mptr) - node(ndilo,mptr) + 1
          ikeep = nx/intratx(level-1)
          lenbc = 2*(ikeep)
c         twice perimeter since saving plus or minus fluxes 
c         plus coarse solution storage
          call reclam(node(ffluxptr,mptr), 2*nvar*lenbc+naux*lenbc)
          mptr  = node(levelptr,mptr)
          if (mptr .ne. 0) go to 40
c
 99    return
       end
