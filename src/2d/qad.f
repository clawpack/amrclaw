c
!> For each coarse-fine interface, a Riemann problem between an inner
!! ghost cell value on the fine grid and cell value in the adjacent coarse
!! cell must be solved and added to corresponding location in
!! **node(ffluxptr, mptr)** for conservative fix later
!!
c -------------------------------------------------------------
c
       subroutine qad(valbig,mitot,mjtot,nvar,
     .                svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy,
     .                maux,aux,auxc1d,delt,mptr)

      use amr_module
       implicit double precision (a-h, o-z)


       logical qprint

       dimension valbig(nvar,mitot,mjtot)
       dimension qc1d(nvar,lenbc)
       dimension svdflx(nvar,lenbc)
       dimension aux(maux,mitot,mjtot)
       dimension auxc1d(maux,lenbc)

c
c ::::::::::::::::::::::::::: QAD ::::::::::::::::::::::::::::::::::
c  are added in to coarse grid value, as a conservation fixup. 
c  Done each fine grid time step. If source terms are present, the
c  coarse grid value is advanced by source terms each fine time step too.

c  No change needed in this sub. for spherical mapping: correctly
c  mapped vals already in bcs on this fine grid and coarse saved
c  vals also properly prepared
c
c Side 1 is the left side of the fine grid patch.  Then go around clockwise.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c      # local storage
c      # note that dimension here are bigger than dimensions used
c      # in rp2, but shouldn't matter since wave is not used in qad
c      # and for other arrays it is only the last parameter that is wrong
c      #  ok as long as meqn, mwaves < maxvar

       integer max1dp1
       dimension ql(nvar,max1d+1),    qr(nvar,max1d+1)
       dimension wave(nvar,mwaves,max1d+1), s(mwaves,max1d+1)
       dimension amdq(nvar,max1d+1),  apdq(nvar,max1d+1)
       dimension auxl(maxaux*(max1d+1)),  auxr(maxaux*(max1d+1))
c
c  WARNING: auxl,auxr dimensioned at max possible, but used as if
c  they were dimensioned as the real maux by max1d+1. Would be better
c  of course to dimension by maux by max1d+1 but this wont work if maux=0
c  So need to access using your own indexing into auxl,auxr.
       iaddaux(iaux,i) = iaux + maux*(i-1)

       data qprint/.false./
c
c      aux is auxiliary array with user parameters needed in Riemann solvers
c          on fine grid corresponding to valbig
c      auxc1d is coarse grid stuff from around boundary, same format as qc1d
c      auxl, auxr are work arrays needed to pass stuff to rpn2
c      maux is the number of aux variables, which may be zero.
c

       max1dp1 = max1d+1

       tgrid = rnode(timemult, mptr)
       if (qprint) 
     .      write(dbugunit,*)" working on grid ",mptr," time ",tgrid
       nc = mjtot-2*nghost
       nr = mitot-2*nghost
       level = node(nestlevel, mptr)
       index = 0

c
c--------
c  side 1
c--------
c
       do j = nghost+1, mjtot-nghost
        if (maux.gt.0) then
          do ma = 1,maux
             if (auxtype(ma).eq."xleft") then
c                # Assuming velocity at left-face, this fix
c                # preserves conservation in incompressible flow:
                 auxl(iaddaux(ma,j-nghost+1)) = aux(ma,nghost+1,j)
               else
c                # Normal case -- we set the aux arrays 
c                # from the cell corresponding  to q
                 auxl(iaddaux(ma,j-nghost+1)) = aux(ma,nghost,j)
               endif
          end do
        endif
        do ivar = 1, nvar
          ql(ivar,j-nghost+1) = valbig(ivar,nghost,j)
        end do
       end do

       lind = 0
       ncrse = (mjtot-2*nghost)/lratioy
       do jc = 1, ncrse
         index = index + 1
         do l = 1, lratioy
         lind = lind + 1
         if (maux.gt.0) then
            do ma=1,maux
               auxr(iaddaux(ma,lind)) = auxc1d(ma,index)
            end do
         endif
         do ivar = 1, nvar
            qr(ivar,lind) = qc1d(ivar,index)
         end do
        end do
       end do
    
       if (qprint) then
         write(dbugunit,*) 'side 1, ql and qr:'
         do i=2,nc
            write(dbugunit,4101) i,qr(1,i-1),ql(1,i)
          enddo
 4101      format(i3,4e16.6)
         if (maux .gt. 0) then
             write(dbugunit,*) 'side 1, auxr:'
             do i=2,nc
                write(dbugunit,4101) i,(auxr(iaddaux(ma,i-1)),ma=1,maux)
                enddo
             write(dbugunit,*) 'side 1, auxl:'
             do i=2,nc
                write(dbugunit,4101) i,(auxl(iaddaux(ma,i)),ma=1,maux)
                enddo
         endif
       endif
 
       call rpn2(1,max1dp1-2*nghost,nvar,mwaves,maux,nghost,
     .              nc+1-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)
c
c we have the wave. for side 1 add into sdflxm
c
       influx = 0
       do j = 1, nc/lratioy
          influx  = influx + 1
          jfine = (j-1)*lratioy
          do ivar = 1, nvar
            do l = 1, lratioy
              svdflx(ivar,influx) = svdflx(ivar,influx) 
     .                     + amdq(ivar,jfine+l+1) * hy * delt
     .                     + apdq(ivar,jfine+l+1) * hy * delt
            end do
          end do
       end do

c--------
c  side 2
c--------
c
       if (mjtot .eq. 2*nghost+1) then
c          # a single row of interior cells only happens when using the
c          # 2d amrclaw code to do a 1d problem with refinement.
c          # (feature added in Version 4.3)
c          # skip over sides 2 and 4 in this case
           go to 299
           endif

       do i = nghost+1, mitot-nghost
        if (maux.gt.0) then
          do ma = 1,maux
             auxr(iaddaux(ma,i-nghost)) = aux(ma,i,mjtot-nghost+1)
          end do 
        endif
        do ivar = 1, nvar
            qr(ivar,i-nghost) = valbig(ivar,i,mjtot-nghost+1)
        end do
       end do  

       lind = 0
       ncrse = (mitot-2*nghost)/lratiox
       do ic = 1, ncrse
         index = index + 1
         do l = 1, lratiox
          lind = lind + 1
          if (maux.gt.0) then
            do ma=1,maux
             if (auxtype(ma).eq."yleft") then
c                # Assuming velocity at bottom-face, this fix
c                # preserves conservation in incompressible flow:
                 ifine = (ic-1)*lratiox + nghost + l
                 auxl(iaddaux(ma,lind+1)) = aux(ma,ifine,mjtot-nghost+1)
              else
                 auxl(iaddaux(ma,lind+1)) = auxc1d(ma,index)
              endif
            end do
          endif
          do ivar = 1, nvar
             ql(ivar,lind+1) = qc1d(ivar,index)
          end do
         end do
       end do
    
       if (qprint) then
         write(dbugunit,*) 'side 2, ql and qr:'
         do i=1,nr
            write(dbugunit,4101) i,ql(1,i+1),qr(1,i)
            enddo
         if (maux .gt. 0) then
             write(dbugunit,*) 'side 2, auxr:'
             do i = 1, nr
                write(dbugunit,4101) i, (auxr(iaddaux(ma,i)),ma=1,maux)
                enddo
             write(dbugunit,*) 'side 2, auxl:'
             do i = 1, nr
                write(dbugunit,4101) i, (auxl(iaddaux(ma,i)),ma=1,maux)
                enddo
         endif
       endif
       call rpn2(2,max1dp1-2*nghost,nvar,mwaves,maux,nghost,
     .              nr+1-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)
c
c we have the wave. for side 2. add into sdflxp
c
       do i = 1, nr/lratiox
          influx  = influx + 1
          ifine = (i-1)*lratiox
          do ivar = 1, nvar
            do l = 1, lratiox
              svdflx(ivar,influx) = svdflx(ivar,influx) 
     .                     - amdq(ivar,ifine+l+1) * hx * delt
     .                     - apdq(ivar,ifine+l+1) * hx * delt
            end do
          end do 
       end do

 299  continue

c--------
c  side 3
c--------
c
       do j = nghost+1, mjtot-nghost
        if (maux.gt.0) then
          do ma = 1,maux
             auxr(iaddaux(ma,j-nghost)) = aux(ma,mitot-nghost+1,j)
          end do
        endif
        do ivar = 1, nvar
           qr(ivar,j-nghost) = valbig(ivar,mitot-nghost+1,j)
        end do
       end do 

       lind = 0
       ncrse = (mjtot-2*nghost)/lratioy
       do jc = 1, ncrse
         index = index + 1
         do l = 1, lratioy
          lind = lind + 1
          if (maux.gt.0) then
            do ma=1,maux
               if (auxtype(ma).eq."xleft") then
c                # Assuming velocity at left-face, this fix
c                # preserves conservation in incompressible flow:
                 jfine = (jc-1)*lratioy + nghost + l
                 auxl(iaddaux(ma,lind+1)) = aux(ma,mitot-nghost+1,jfine)
               else
                 auxl(iaddaux(ma,lind+1)) = auxc1d(ma,index)
               endif
             end do
          endif
          do ivar = 1, nvar
              ql(ivar,lind+1) = qc1d(ivar,index)
          end do
         end do
       end do
    
       if (qprint) then
         write(dbugunit,*) 'side 3, ql and qr:'
         do i=1,nc
            write(dbugunit,4101) i,ql(1,i),qr(1,i)
            enddo
       endif
       call rpn2(1,max1dp1-2*nghost,nvar,mwaves,maux,nghost,
     .              nc+1-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)
c
c we have the wave. for side 3 add into sdflxp
C
       do j = 1, nc/lratioy
          influx  = influx + 1
          jfine = (j-1)*lratioy
          do ivar = 1, nvar
            do l = 1, lratioy
              svdflx(ivar,influx) = svdflx(ivar,influx) 
     .                     - amdq(ivar,jfine+l+1) * hy * delt
     .                     - apdq(ivar,jfine+l+1) * hy * delt
             end do
           end do
        end do

c--------
c  side 4
c--------
c
       if (mjtot .eq. 2*nghost+1) then
c          # a single row of interior cells only happens when using the
c          # 2d amrclaw code to do a 1d problem with refinement.
c          # (feature added in Version 4.3)
c          # skip over sides 2 and 4 in this case
           go to 499
           endif
c
       do i = nghost+1, mitot-nghost
        if (maux.gt.0) then
          do ma = 1,maux
             if (auxtype(ma).eq."yleft") then
c                # Assuming velocity at bottom-face, this fix
c                # preserves conservation in incompressible flow:
                 auxl(iaddaux(ma,i-nghost+1)) = aux(ma,i,nghost+1)
              else
                 auxl(iaddaux(ma,i-nghost+1)) = aux(ma,i,nghost)
              endif
           end do
        endif
        do ivar = 1, nvar
          ql(ivar,i-nghost+1) = valbig(ivar,i,nghost)
        end do
       end do

       lind = 0
       ncrse = (mitot-2*nghost)/lratiox
       do ic = 1, ncrse
         index = index + 1
         do l = 1, lratiox
          lind = lind + 1
          if (maux.gt.0) then
             do ma=1,maux
                auxr(iaddaux(ma,lind)) = auxc1d(ma,index)
             end do
          endif
          do ivar = 1, nvar
              qr(ivar,lind) = qc1d(ivar,index)
          end do
         end do
       end do
    
       if (qprint) then
         write(dbugunit,*) 'side 4, ql and qr:'
         do i=1,nr
            write(dbugunit,4101) i, ql(1,i),qr(1,i)
            enddo
       endif
       call rpn2(2,max1dp1-2*nghost,nvar,mwaves,maux,nghost,
     .              nr+1-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)
c
c we have the wave. for side 4. add into sdflxm
c
       do i = 1, nr/lratiox
          influx  = influx + 1
          ifine = (i-1)*lratiox
          do ivar = 1, nvar
            do l = 1, lratiox
              svdflx(ivar,influx) = svdflx(ivar,influx) 
     .                     + amdq(ivar,ifine+l+1) * hx * delt
     .                     + apdq(ivar,ifine+l+1) * hx * delt
             end do
           end do
        end do

 499   continue

c      # for source terms:
       if (method(5) .ne. 0) then   ! should I test here if index=0 and all skipped?
           call src1d(nvar,nghost,lenbc,qc1d,maux,auxc1d,tgrid,delt)
c      # how can this be right - where is the integrated src term used?
           endif
           
       return
       end
