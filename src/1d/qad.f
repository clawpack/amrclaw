c
c -------------------------------------------------------------
c
       subroutine qad(valbig,mitot,nvar,
     .                svdflx,qc1d,lenbc,lratiox,hx,
     .                maux,aux,auxc1d,delt,mptr)

      use amr_module
       implicit double precision (a-h, o-z)


       logical qprint

       dimension valbig(nvar,mitot)
       dimension qc1d(nvar,lenbc)
       dimension svdflx(nvar,lenbc)
       dimension aux(maux,mitot)
       dimension auxc1d(maux,lenbc)

c
c ::::::::::::::::::::::::::: QAD ::::::::::::::::::::::::::::::::::
c  solve RP between ghost cell value on fine grid and coarse grid
c  value that ghost cell overlaps. The resulting fluctuations
c  are added in to coarse grid value, as a conservation fixup. 
c  Done each fine grid time step. If source terms are present, the
c  coarse grid value is advanced by source terms each fine time step too.

c  No change needed in this sub. for spherical mapping: correctly
c  mapped vals already in bcs on this fine grid and coarse saved
c  vals also properly prepared
c
c Side 1 is the left side of the fine grid patch.
c Side 2 is the right side of the fine grid patch.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c      # local storage
c      # note that dimension here are bigger than dimensions used
c      # in rp1, but shouldn't matter since wave is not used in qad
c      # and for other arrays it is only the last parameter that is wrong
c      #  ok as long as meqn, mwaves < maxvar

       integer  max1dp1
       dimension ql(nvar,max1d+1),    qr(nvar,max1d+1)
       dimension wave(nvar,mwaves,max1d+1), s(mwaves,max1d+1)
       dimension amdq(nvar,max1d+1),  apdq(nvar,max1d+1)
       dimension auxl(maxaux*(max1d+1)),  auxr(maxaux*(max1d+1))
c
c  WARNING: auxl,auxr dimensioned at max possible, but used as if
c  they were dimensioned as the real maux by max1dp1. Would be better
c  of course to dimension by maux by max1dp1 but this wont work if maux=0
c  So need to access using your own indexing into auxl,auxr.
       iaddaux(iaux,i) = iaux + maux*(i-1)

       data qprint/.false./
c
c      aux is auxiliary array with user parameters needed in Riemann solvers
c          on fine grid corresponding to valbig
c      auxc1d is coarse grid stuff from around boundary, same format as qc1d
c      auxl, auxr are work arrays needed to pass stuff to rp1
c      maux is the number of aux variables, which may be zero.
c
       max1dp1 = max1d + 1

       tgrid = rnode(timemult, mptr)
       if (qprint) 
     .      write(dbugunit,*)" working on grid ",mptr," time ",tgrid
       level = node(nestlevel, mptr)
       index = 0

c
c--------
c  side 1
c--------
c
       if (maux.gt.0) then
          do 5 ma = 1,maux
             if (auxtype(ma).eq."xleft") then
c                # Assuming velocity at left-face, this fix
c                # preserves conservation in incompressible flow:
                 auxl(iaddaux(ma,2)) = aux(ma,nghost+1)
               else
c                # Normal case -- we set the aux arrays 
c                # from the cell corresponding  to q
                 auxl(iaddaux(ma,2)) = aux(ma,nghost)
               endif
  5          continue
          endif
       do 10 ivar = 1, nvar
         ql(ivar,2) = valbig(ivar,nghost)
 10    continue

       index = index + 1
       lind = 1
       if (maux.gt.0) then
         do 24 ma=1,maux
            auxr(iaddaux(ma,lind)) = auxc1d(ma,index)
   24       continue
       endif
       do 25 ivar = 1, nvar
         qr(ivar,lind) = qc1d(ivar,index)
 25    continue
    
       if (qprint) then
         write(dbugunit,*) 'side 1, ql and qr:'
            write(dbugunit,4101) i,qr(1,1),ql(1,2)
 4101      format(i3,4e16.6)
         if (maux .gt. 0) then
             write(dbugunit,*) 'side 1, auxr:'
             write(dbugunit,4101) i,(auxr(iaddaux(ma,1)),ma=1,maux)
             write(dbugunit,*) 'side 1, auxl:'
             write(dbugunit,4101) i,(auxl(iaddaux(ma,2)),ma=1,maux)
         endif
       endif

       call rp1(max1dp1-2*nghost,nvar,mwaves,maux,nghost,
     .              2-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)


c
c we have the wave. for side 1 add into sdflxm
c
       influx  = 1
       do 40 ivar = 1, nvar
           svdflx(ivar,influx) = svdflx(ivar,influx)
     .                     + amdq(ivar,2) * delt
     .                     + apdq(ivar,2) * delt
 40    continue

c--------
c  side 2
c--------
c
        if (maux.gt.0) then
          do 205 ma = 1,maux
             auxr(iaddaux(ma,1)) = aux(ma,mitot-nghost+1)
 205         continue
          endif
        do 210 ivar = 1, nvar
          qr(ivar,1) = valbig(ivar,mitot-nghost+1)
 210      continue

         index = index + 1
         lind = 1
         if (maux.gt.0) then
            do 224 ma=1,maux
             if (auxtype(ma).eq."xleft") then
c                # Assuming velocity at left-face, this fix
c                # preserves conservation in incompressible flow:
                 auxl(iaddaux(ma,lind+1)) = aux(ma,mitot-nghost+1)
               else
                 auxl(iaddaux(ma,lind+1)) = auxc1d(ma,index)
               endif
  224          continue
            endif
         do 225 ivar = 1, nvar
             ql(ivar,lind+1) = qc1d(ivar,index)
 225     continue

       if (qprint) then
         write(dbugunit,*) 'side 2, ql and qr:'
            write(dbugunit,4101) i,ql(1,2),qr(1,1)
       endif

       call rp1(max1dp1-2*nghost,nvar,mwaves,maux,nghost,
     .              2-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)

c
c we have the wave. for side 2 add into sdflxp
C
          influx  = influx + 1
          do 340 ivar = 1, nvar
              svdflx(ivar,influx) = svdflx(ivar,influx)
     .                     - amdq(ivar,2) * delt
     .                     - apdq(ivar,2) * delt
 340       continue

c      # for source terms:
       if (method(4) .ne. 0) then   ! should I test here if index=0 and all skipped?
           call src1d(nvar,nghost,lenbc,qc1d,maux,auxc1d,tgrid,delt)
c      # how can this be right - where is the integrated src term used?
           endif

       return
       end
