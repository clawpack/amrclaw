c
c -----------------------------------------------------------------
c
      subroutine setgrd (nvar,cut,naux,dtinit,start_time)
c
      use amr_module
      implicit double precision (a-h,o-z)

      logical  vtime
      data     vtime/.false./

c     # may as well not bother to calculate time step for error est.
c
c :::::::::::::::::::::::::::: SETGRD :::::::::::::::::::::::::::::::;
c  set up the entire tree/grid structure.  only at this time t = 0
c  can we take advantage of initialization routines.
c  remember that regridding/error estimation needs to have two
c  time steps of soln. values.
c  6/21/05: added dtinit arg. to allow for better choice of initial timestep
c   as discovered by advance/setgrd in first step.
c ::::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::::::::;
c
      dtinit = possk(1)
      if (mxnest .eq. 1) go to 99
c
      levnew =  2
      time   = start_time
      verbosity_regrid = method(4)
c
 10   if (levnew .gt. mxnest) go to 30
          levold = levnew - 1
          if (lstart(levold) .eq. 0) go to 30
          lbase  = levold
          lfnew  = lbase
c
c  set up level to be flagged. need a solution t=0,and if
c  using richardson  error estimation then one at dt too.
c  then richardson makes one at 2*dt.  for just flag2refine 
c  only need to set boundary conditions, done from grdfit->flagger.
c
          if (flag_richardson) then
             call advanc(levold,nvar,dtlev,vtime,naux)
             evol = evol + rvol !time stepping stats for error estimation
             rvol = 0.d0        ! reset since not 'real' time stepping stats
             kfac = 1
             do i = 1, levold-1
                kfac = kfac * kratio(i)
             end do
             dtinit = min(dtinit, dtlev*kfac)
             
c            dont count it in real integration stats
             do 20 level=1,mxnest
 20             rvoll(level) = 0.d0
             endif
c
c  flag, cluster, and make new grids. grdfit set bcs, controls flagging,
c  colating and making grids. But advanc called above since if using
c  richardson, it assumes two levels of solution already exist
c
         call grdfit(lbase,levold,nvar,naux,cut,time,start_time)
         if (newstl(levnew) .ne. 0) lfnew = levnew
c
c  init new level. after each iteration. fix the data structure
c  also reinitalize coarser grids so fine grids can be advanced
c  and interpolate correctly for their bndry vals from coarser grids.
c
         call ginit(newstl(levnew),.true., nvar, naux,start_time)
         lstart(levnew) = newstl(levnew)
         lfine = lfnew
         call ginit(lstart(levold),.false., nvar, naux,start_time)
c
c count number of grids on newly created levels (needed for openmp
c parallelization). this is also  done in regridding.
c  set up numgrids now for new level, since advanc will use it for parallel execution
c 
         mptr = lstart(levnew)
         ngrids = 0
         ncells = 0
         do while (mptr .gt. 0)
            ngrids = ngrids + 1
            ncells = ncells + (node(ndihi,mptr)-node(ndilo,mptr)+1)
     .                      * (node(ndjhi,mptr)-node(ndjlo,mptr)+1)
            mptr = node(levelptr, mptr)
          end do
          numgrids(levnew) = ngrids
          numcells(levnew) = ncells
          avenumgrids(levnew) = avenumgrids(levnew) + ngrids
          iregridcount(levnew) = iregridcount(levnew) + 1
          if (ngrids .gt. 1) call arrangeGrids(levnew,ngrids)
          if (verbosity_regrid .ge. levnew) then
             write(*,100) ngrids,ncells,levnew
 100         format("there are ",i4," grids with ",i8,
     &               " cells at level ", i3)
          endif 
c
c     need to make gridList here before calling again to make finer grids.
c     This is because ths list if used in advanc. level 1 is called from domain
      call makeGridList(lbase)
      call makeBndryList(levnew)
c
      levnew = levnew + 1
      go to 10
 30   continue
c
c  switch location of old and new storage for soln. vals, 
c  and reset time to 0.0 (or initial time start_time)
c
c      if (mxnest .eq. 1) go to 99  shouldnt be nec. tested above.
c
      lev = 1
 40   if ((lev .eq. mxnest) .or. (lev .gt. lfine))  go to 60
        mptr = lstart(lev)
 50        itemp                = node(store1,mptr)
           node(store1,mptr)    = node(store2,mptr)
           node(store2,mptr)    = itemp
           rnode(timemult,mptr) = start_time
           mptr                 = node(levelptr,mptr)
           if (mptr .ne. 0) go to 50
       lev = lev + 1
       go to 40
 60   continue
c
c initial updating so can do conservation check. can do before
c bndry flux arrays set, since dont need them for this
c
      do 65 level = 1, lfine-1
         call update(lfine-level,nvar,naux)
 65   continue
c
c set up boundary flux conservation arrays
c
      do 70 level = 1, lfine-1
         call prepf(level+1,nvar,naux)
         call prepc(level,nvar)
 70   continue
c
c
c
c
 99   continue
      return
      end
