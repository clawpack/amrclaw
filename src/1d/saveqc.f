c
c  ================================================================
      subroutine saveqc(level,nvar,naux)
c  ================================================================
c
      use amr_module
      implicit double precision (a-h,o-z)
      
      !for setaux timing
      integer :: clock_start, clock_finish, clock_rate
      real(kind=8) :: cpu_start, cpu_finish
      

      logical sticksout, found
!     make fliparray largest possible grid size
      dimension fliparray(2*max1d*nghost*(nvar+naux))
c
c ::::::::::::::::::::::::: SAVEQC :::::::::::::::::::::::::::::::::
c  prepare new fine grids to save fluxes after each integration step
c  for future conservative fix-up on coarse grids.
c  save all boundary fluxes of fine grid (even if on a  phys. bndry.) -
c  but only save space for every intrat of them. 
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      levc = level - 1
      hxc  = hxposs(levc)
      ng   = 0  ! no ghost cells on coarsened enlarged patch

      mkid = lstart(level)
 10   if (mkid .eq. 0) go to 99
          nx    = node(ndihi,mkid)-node(ndilo,mkid) + 1
          lenbc = 2
          ist   = node(ffluxptr,mkid)
          time = rnode(timemult,mkid)

c         make coarsened enlarged patch for conservative fixup
          ilo = node(ndilo,mkid)
          ihi = node(ndihi,mkid)
          iclo = ilo/intratx(level-1) - 1
          ichi = (ihi+1)/intratx(level-1)
          nrow = ichi-iclo+1
          xl   = rnode(cornxlo,mkid) - hxc
          xr   = rnode(cornxhi,mkid) + hxc

          loctmp = igetsp(nrow*(nvar+naux))
          loctx  = loctmp + nrow*nvar
          do i = 1, nrow*naux
             alloc(loctx+i-1) = NEEDS_TO_BE_SET
          end do
          locaux = node(storeaux,mkid)

          if (iclo .lt. 0 .or. ichi .eq. iregsz(levc)) then
            sticksout = .true.
          else
            sticksout = .false.
          endif

          if (sticksout .and. xperdom) then
             !iperim = nrow
             !locflip = igetsp(iperim*nghost*(nvar+naux))

             call preicall(alloc(loctmp),alloc(loctx),nrow,nvar,
     .                     naux,iclo,ichi,level-1,
     .                     fliparray)
!     .                     alloc(locflip))
!             call reclam(locflip,iperim*nghost*(nvar+naux))
          else
             call icall(alloc(loctmp),alloc(loctx),nrow,nvar,naux,
     .                   iclo,ichi,level-1,1)
          endif
!         in case any part sticks out of domain still need to set remaining aux
!         cells
          if (naux .gt. 0 .and. sticksout) then  
             call system_clock(clock_start,clock_rate)
             call cpu_time(cpu_start)
             call setaux(ng,nrow,xl,hxc,naux,alloc(loctx))
             call system_clock(clock_finish,clock_rate)
             call cpu_time(cpu_finish)
             timeSetaux = timeSetaux + clock_finish - clock_start
             timeSetauxCPU = timeSetauxCPU + cpu_finish - cpu_start
          endif
!--          found = .false.
!--          do i = 1, naux*nrow, naux
!--             if (alloc(loctx+i-1) .eq. NEEDS_TO_BE_SET) then
!--                 found = .true.
!--             endif
!--          end do
!--          if (found)  write(*,*) "still have unset aux vals in qad"
          call bc1amr(alloc(loctmp),alloc(loctx),nrow,nvar,naux,
     .                hxc,level,time,
     .                xl,xr)
          call cstore(alloc(loctmp),nrow,nvar,
     .                alloc(ist+nvar*lenbc),lenbc,naux,alloc(loctx),
     .                alloc(ist+2*nvar*lenbc))
          call reclam(loctmp,nrow*(nvar+naux))

          mkid = node(levelptr,mkid)
          go to 10
 99    return
       end
