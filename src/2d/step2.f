c
c
c     ==========================================================
      subroutine step2(maxm,maxmx,maxmy,meqn,maux,mbc,mx,my,
     &                 qold,aux,dx,dy,dt,cflgrid,
     &                 fm,fp,gm,gp,rpn2,rpt2)
c     ==========================================================
c
c     # clawpack routine ...  modified for AMRCLAW
c
c     # Take one time step, updating q.
c     # On entry, qold gives
c     #    initial data for this step
c     #    and is unchanged in this version.
c    
c     # fm, fp are fluxes to left and right of single cell edge
c     # See the flux2 documentation for more information.
c
c
      use amr_module
      implicit double precision (a-h,o-z)
      external rpn2, rpt2

      dimension qold(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension aux(maux,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension   fm(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension   fp(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension   gm(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension   gp(meqn,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)

      ! These are stack based here such that each thread can have their own copy
      dimension faddm(meqn,1-mbc:maxm+mbc)
      dimension faddp(meqn,1-mbc:maxm+mbc)
      dimension gaddm(meqn,1-mbc:maxm+mbc,2)
      dimension gaddp(meqn,1-mbc:maxm+mbc,2)
      dimension  q1d(meqn,1-mbc:maxm+mbc)
      dimension aux1(maux,1-mbc:maxm+mbc)
      dimension aux2(maux,1-mbc:maxm+mbc)
      dimension aux3(maux,1-mbc:maxm+mbc)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      
      dimension  wave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension     s(mwaves, 1-mbc:maxm + mbc)
      dimension  amdq(meqn,1-mbc:maxm + mbc)
      dimension  apdq(meqn,1-mbc:maxm + mbc)
      dimension  cqxx(meqn,1-mbc:maxm + mbc)
      dimension bmadq(meqn,1-mbc:maxm + mbc)
      dimension bpadq(meqn,1-mbc:maxm + mbc)
      
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
c     # store mesh parameters that may be needed in Riemann solver but not
c     # passed in...
      dxcom = dx
      dycom = dy
      dtcom = dt
c
c
c     # partition work array into pieces needed for local storage in
c     # flux2 routine.  Find starting index of each piece:
c
C       i0wave = 1
C       i0s = i0wave + (maxm+2*mbc)*meqn*mwaves
C       i0amdq = i0s + (maxm+2*mbc)*mwaves
C       i0apdq = i0amdq + (maxm+2*mbc)*meqn
C       i0cqxx = i0apdq + (maxm+2*mbc)*meqn
C       i0bmadq = i0cqxx + (maxm+2*mbc)*meqn
C       i0bpadq = i0bmadq + (maxm+2*mbc)*meqn
C       iused = i0bpadq + (maxm+2*mbc)*meqn - 1
c
C       if (iused.gt.mwork) then
C c        # This shouldn't happen due to checks in claw2
C          write(outunit,*) 'not enough work space in step2'
C          write(*      ,*) 'not enough work space in step2'
C          stop 
C          endif
c
c
      cflgrid = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
c
      fm = 0.d0
      fp = 0.d0
      gm = 0.d0
      gp = 0.d0
c
      if (mcapa == 0) then
c         No capa array:
          dtdx1d = dtdx
          dtdy1d = dtdy
      endif
c
c
c     # perform x-sweeps
c     ==================
c
#ifdef SWEEP_THREADING
!$OMP PARALLEL DO PRIVATE(j,i,m,ma,dtdx,jcom)
!$OMP&            PRIVATE(faddm,faddp,gaddm,gaddp,q1d,dtdx1d)
!$OMP&            PRIVATE(aux1,aux2,aux3)
!$OMP&            PRIVATE(wave,s,amdq,apdq,cqxx,bmadq,bpadq)
!$OMP&            PRIVATE(cfl1d)
!$OMP&            SHARED(mx,my,maxm,maux,mcapa,mbc,meqn)
!$OMP&            SHARED(cflgrid,fm,fp,gm,gp,qold,aux)
!$OMP&            DEFAULT(none)
#endif
      do 50 j = 0,my+1
         if (my.eq.1 .and. j.ne.1) go to 50  !# for 1d AMR
c
c        # copy data along a slice into 1d arrays:
         do 20 i = 1-mbc, mx+mbc
           do 20 m=1,meqn
               q1d(m,i) = qold(m,i,j)
   20          continue
c
         if (mcapa.gt.0)  then
           do 21 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(mcapa,i,j)
   21          continue
           endif
c
         if (maux .gt. 0)  then
            do 22 i = 1-mbc, mx+mbc
              do 22 ma=1,maux
                 aux1(ma,i) = aux(ma,i,j-1)
                 aux2(ma,i) = aux(ma,i,j  )
                 aux3(ma,i) = aux(ma,i,j+1)
   22            continue
           endif
c
c
c        # Store the value of j along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         jcom = j  
c                  
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(1,maxm,meqn,maux,mbc,mx,
     &              q1d,dtdx1d,aux1,aux2,aux3,
     &              faddm,faddp,gaddm,gaddp,cfl1d,
     &              wave,s,amdq,apdq,
     &              cqxx,bmadq,bpadq,rpn2,rpt2)
#ifdef SWEEP_THREADING
!$OMP CRITICAL (cfl_row)
#endif
         cflgrid = dmax1(cflgrid,cfl1d)
#ifdef SWEEP_THREADING
!$OMP END CRITICAL (cfl_row)
#endif
c
c        # update fluxes for use in AMR:
c
#ifdef SWEEP_THREADING
!$OMP CRITICAL (flux_accumulation)
#endif
         do m=1,meqn
            do i=1,mx+1
               fm(m,i,j) = fm(m,i,j) + faddm(m,i)
               fp(m,i,j) = fp(m,i,j) + faddp(m,i)
               gm(m,i,j) = gm(m,i,j) + gaddm(m,i,1)
               gp(m,i,j) = gp(m,i,j) + gaddp(m,i,1)
               gm(m,i,j+1) = gm(m,i,j+1) + gaddm(m,i,2)
               gp(m,i,j+1) = gp(m,i,j+1) + gaddp(m,i,2)
            enddo
         enddo
#ifdef SWEEP_THREADING
!$OMP END CRITICAL (flux_accumulation)
#endif
   50    continue
#ifdef SWEEP_THREADING
!$OMP END PARALLEL DO
#endif
c
c
c
c     # perform y sweeps
c     ==================
c
c

      if (my.eq.1) go to 101   !# for 1d AMR
#ifdef SWEEP_THREADING
!$OMP PARALLEL DO PRIVATE(j,i,m,ma,dtdy,icom)
!$OMP&            PRIVATE(faddm,faddp,gaddm,gaddp,q1d,dtdy1d)
!$OMP&            PRIVATE(aux1,aux2,aux3)
!$OMP&            PRIVATE(wave,s,amdq,apdq,cqxx,bmadq,bpadq)
!$OMP&            PRIVATE(cfl1d)
!$OMP&            SHARED(mx,my,maxm,maux,mcapa,mbc,meqn)
!$OMP&            SHARED(cflgrid,fm,fp,gm,gp,qold,aux)
!$OMP&            DEFAULT(none)
#endif
      do 100 i = 0, mx+1
c
c        # copy data along a slice into 1d arrays:
         do 70 j = 1-mbc, my+mbc
           do 70 m=1,meqn
               q1d(m,j) = qold(m,i,j)
   70          continue
c
         if (mcapa.gt.0)  then
           do 71 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(mcapa,i,j)
   71          continue
           endif
c
         if (maux .gt. 0)  then
            do 72 j = 1-mbc, my+mbc
              do 72 ma=1,maux
                 aux1(ma,j) = aux(ma,i-1,j)
                 aux2(ma,j) = aux(ma,i,  j)
                 aux3(ma,j) = aux(ma,i+1,j)
   72            continue
           endif
c
c
c        # Store the value of i along this slice in the common block
c        # comxyt in case it is needed in the Riemann solver (for
c        # variable coefficient problems)
         icom = i  
c                  
c        # compute modifications fadd and gadd to fluxes along this slice:
         call flux2(2,maxm,meqn,maux,mbc,my,
     &              q1d,dtdy1d,aux1,aux2,aux3,
     &              faddm,faddp,gaddm,gaddp,cfl1d,
     &              wave,s,amdq,apdq,
     &              cqxx,bmadq,bpadq,rpn2,rpt2)
c

#ifdef SWEEP_THREADING
!$OMP CRITICAL (cfl_row)
#endif
           cflgrid = dmax1(cflgrid,cfl1d)
#ifdef SWEEP_THREADING
!$OMP END CRITICAL (cfl_row)
#endif

c
c        # 
c        # update fluxes for use in AMR:
c
#ifdef SWEEP_THREADING
!$OMP CRITICAL (flux_accumulation)
#endif
         do j=1,my+1
            do m=1,meqn
               gm(m,i,j) = gm(m,i,j) + faddm(m,j)
               gp(m,i,j) = gp(m,i,j) + faddp(m,j)
               fm(m,i,j) = fm(m,i,j) + gaddm(m,j,1)
               fp(m,i,j) = fp(m,i,j) + gaddp(m,j,1)
               fm(m,i+1,j) = fm(m,i+1,j) + gaddm(m,j,2)
               fp(m,i+1,j) = fp(m,i+1,j) + gaddp(m,j,2)
            enddo
         enddo
#ifdef SWEEP_THREADING
!$OMP END CRITICAL (flux_accumulation)
#endif
  100    continue
#ifdef SWEEP_THREADING
!$OMP END PARALLEL DO
#endif
c
  101 continue
c
      return
      end
