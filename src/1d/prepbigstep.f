c
c-------------------------------------------------------------------------------------
c
       subroutine prepbigstep(nvar,naux,lcheck,mptr,nx,midub,
     .                       valbgc,auxbgc,mi2tot)

       use amr_module
       implicit double precision (a-h,o-z)

       double precision valdub(nvar,midub)
       double precision auxdub(naux,midub)
       double precision valbgc(nvar,mi2tot)
       double precision auxbgc(naux,mi2tot)
       dimension fp(nvar,mi2tot)
       dimension fm(nvar,mi2tot)
       
       !for setaux timing
       integer :: clock_start, clock_finish, clock_rate
       real(kind=8) :: cpu_start, cpu_finish

          hx  = hxposs(lcheck)
          hx2 = 2.d0*hx
          dt  = possk(lcheck)
          dt2 = 2. * dt
          time  = rnode(timemult,mptr)
          tpre  = time - dt

          mitot  = nx + 2*nghost
          ng2    = 2*nghost
          locold = node(store2,mptr)
          xlow   = rnode(cornxlo,mptr) - nghost*hx2
c

c         # transfer soln. into grid with twice the ghost cells
          call copysol(valdub,alloc(locold),nvar,mitot,
     1              nghost,midub,ng2)

c
          if (naux .gt. 0) then
              xl     = rnode(cornxlo, mptr)
              mx = midub - 4*nghost
              auxdub = NEEDS_TO_BE_SET  ! signal that needs a val
              
              call system_clock(clock_start, clock_rate)
              call cpu_time(cpu_start)
              call setaux(2*nghost,mx,xl,hx,
     &                    naux,auxdub)
              call system_clock(clock_finish, clock_rate)
              call cpu_time(cpu_finish)
              timeSetaux = timeSetaux + clock_finish - clock_start
              timeSetauxCPU = timeSetauxCPU + cpu_finish - cpu_start
              
              call auxcoarsen(auxdub,midub,
     1                     auxbgc,mi2tot,naux,auxtype)
          endif

c         # fill it - use enlarged (before coarsening) aux arrays
          call bound(tpre,nvar,ng2,valdub,midub,mptr,
     1               auxdub,naux)

c         coarsen by 2 in every direction
          call coarsen(valdub,midub,
     1                 valbgc,mi2tot,nvar)

          call stepgrid(valbgc,fm,fp,
     1                mi2tot,nghost,
     2                dt2,dtnew2,hx2,nvar,
     3                xlow,tpre,mptr,naux,auxbgc)

c         update counts for error estimation work
          evol = evol + (nx/2)
           return
           end
