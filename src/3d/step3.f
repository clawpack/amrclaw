c
c     ==================================================================
      subroutine step3(maxm,maxmx,maxmy,maxmz,meqn,maux,mbc,mx,my,mz,
     &                qold,aux,dx,dy,dz,dt,cflgrid,
     &                 fm,fp,gm,gp,hm,hp,
     &                 faddm,faddp,gadd,hadd,
     &                 q1d,dtdx1d,dtdy1d,dtdz1d,
     &                 aux1,aux2,aux3,work,mwork,rpn3,rpt3,rptt3)
c     ==================================================================

c     # clawpack routine ...  modified for AMRCLAW

c
c     # Take one time step, updating q.
c     # On entry, qold gives
c     #    initial data for this step
c     #    and is unchanged in this version.
c
c     # fm, fp are fluxes to left and right of single cell edge
c
c     # See the flux3 documentation for more information.
c
c
      use amr_module
      implicit double precision(a-h,o-z)

      external rpn3,rpt3,rptt3

      dimension
     & qold(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 1-mbc:maxmz+mbc)
      dimension q1d(meqn, 1-mbc:maxm+mbc)
c
      dimension fm(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)
      dimension fp(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)
      dimension gm(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)
      dimension gp(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)
      dimension hm(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)
      dimension hp(meqn,1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,1-mbc:maxmz+mbc)
c
      dimension faddm(meqn, 1-mbc:maxm+mbc)
      dimension faddp(meqn, 1-mbc:maxm+mbc)
      dimension  gadd(meqn, 1-mbc:maxm+mbc, 2, -1:1)
      dimension  hadd(meqn, 1-mbc:maxm+mbc, 2, -1:1)
c
      dimension
     & aux(maux, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 1-mbc:maxmz+mbc)
      dimension aux1(maux, 1-mbc:maxm+mbc, 3)
      dimension aux2(maux, 1-mbc:maxm+mbc, 3)
      dimension aux3(maux, 1-mbc:maxm+mbc, 3)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      dimension dtdz1d(1-mbc:maxm+mbc)
      dimension work(mwork)
c$$$      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom
c
c     # store mesh parameters that may be needed in Riemann solver but not
c     # passed in...
c$$$      dxcom = dx
c$$$      dycom = dy
c$$$      dtcom = dt
c
c
c     # partition work array into pieces needed for local storage in
c     # flux3 routine.  Find starting index of each piece:
c
      i0wave     = 1
      i0s        = i0wave     + (maxm+2*mbc)*meqn*mwaves
      i0amdq     = i0s        + (maxm+2*mbc)*mwaves
      i0apdq     = i0amdq     + (maxm+2*mbc)*meqn
      i0cqxx     = i0apdq     + (maxm+2*mbc)*meqn
      i0bmamdq   = i0cqxx     + (maxm+2*mbc)*meqn
      i0bmapdq   = i0bmamdq   + (maxm+2*mbc)*meqn
      i0bpamdq   = i0bmapdq   + (maxm+2*mbc)*meqn
      i0bpapdq   = i0bpamdq   + (maxm+2*mbc)*meqn
      i0cmamdq   = i0bpapdq   + (maxm+2*mbc)*meqn
      i0cmapdq   = i0cmamdq   + (maxm+2*mbc)*meqn
      i0cpamdq   = i0cmapdq   + (maxm+2*mbc)*meqn
      i0cpapdq   = i0cpamdq   + (maxm+2*mbc)*meqn
      i0cmamdq2  = i0cpapdq   + (maxm+2*mbc)*meqn
      i0cmapdq2  = i0cmamdq2  + (maxm+2*mbc)*meqn
      i0cpamdq2  = i0cmapdq2  + (maxm+2*mbc)*meqn
      i0cpapdq2  = i0cpamdq2  + (maxm+2*mbc)*meqn

      i0bmcqxxm   = i0cpapdq2  + (maxm+2*mbc)*meqn
      i0bmcqxxp   = i0bmcqxxm  + (maxm+2*mbc)*meqn
      i0bpcqxxm   = i0bmcqxxp  + (maxm+2*mbc)*meqn
      i0bpcqxxp   = i0bpcqxxm  + (maxm+2*mbc)*meqn
      i0cmcqxxm   = i0bpcqxxp  + (maxm+2*mbc)*meqn
      i0cmcqxxp   = i0cmcqxxm  + (maxm+2*mbc)*meqn
      i0cpcqxxm   = i0cmcqxxp  + (maxm+2*mbc)*meqn
      i0cpcqxxp   = i0cpcqxxm  + (maxm+2*mbc)*meqn

      work(i0bmcqxxm) = i0bmcqxxm
      work(i0bmcqxxp) = i0bmcqxxp
      work(i0bpcqxxm) = i0bpcqxxm
      work(i0bpcqxxp) = i0bpcqxxp
      work(i0cmcqxxm) = i0cmcqxxm
      work(i0cmcqxxp) = i0cmcqxxp
      work(i0cpcqxxm) = i0cpcqxxm
      work(i0cpcqxxp) = i0cpcqxxp

      i0bmcmamdq = i0cpcqxxp  + (maxm+2*mbc)*meqn
      i0bmcmapdq = i0bmcmamdq + (maxm+2*mbc)*meqn
      i0bpcmamdq = i0bmcmapdq + (maxm+2*mbc)*meqn
      i0bpcmapdq = i0bpcmamdq + (maxm+2*mbc)*meqn
      i0bmcpamdq = i0bpcmapdq + (maxm+2*mbc)*meqn
      i0bmcpapdq = i0bmcpamdq + (maxm+2*mbc)*meqn
      i0bpcpamdq = i0bmcpapdq + (maxm+2*mbc)*meqn
      i0bpcpapdq = i0bpcpamdq + (maxm+2*mbc)*meqn
      iused      = i0bpcpapdq + (maxm+2*mbc)*meqn - 1
c
      if (iused.gt.mwork) then
c        # This shouldnt happen if mwork is set properly in stepgrid3
         write(6,*) '*** not enough work space in step3'
         write(6,*) '*** check parameter mwork in stepgrid3'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
         stop
      endif
c
c
      cflgrid = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
      dtdz = dt/dz
c
      do 10 i=1-mbc,mx+mbc
         do 10 j=1-mbc,my+mbc
            do 10 k=1-mbc,mz+mbc
               do 10 m=1,meqn
                  fm(m,i,j,k) = 0.d0
                  fp(m,i,j,k) = 0.d0
                  gm(m,i,j,k) = 0.d0
                  gp(m,i,j,k) = 0.d0
                  hm(m,i,j,k) = 0.d0
                  hp(m,i,j,k) = 0.d0
 10            continue

c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
            dtdz1d(i) = dtdz
    5       continue
         endif
c
c
c     # perform x-sweeps
c     ==================
c
      do 50 k = 0,mz+1
         do 50 j = 0,my+1
c
            do 20 i = 1-mbc, mx+mbc
               do 20 m=1,meqn
c                 # copy data along a slice into 1d array:
                  q1d(m,i) = qold(m,i,j,k)
   20          continue
c
         if (mcapa.gt.0)  then
           do 23 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(mcapa,i,j,k)
   23      continue
         endif
c
         if (maux .gt. 0)  then
            do 22 i = 1-mbc, mx+mbc
               do 22 ma=1,maux
                 aux1(ma,i,1) = aux(ma,i,j-1,k-1)
                 aux1(ma,i,2) = aux(ma,i,j-1,k)
                 aux1(ma,i,3) = aux(ma,i,j-1,k+1)
                 aux2(ma,i,1) = aux(ma,i,j,k-1)
                 aux2(ma,i,2) = aux(ma,i,j,k)
                 aux2(ma,i,3) = aux(ma,i,j,k+1)
                 aux3(ma,i,1) = aux(ma,i,j+1,k-1)
                 aux3(ma,i,2) = aux(ma,i,j+1,k)
                 aux3(ma,i,3) = aux(ma,i,j+1,k+1)
   22          continue
           endif
c
c           # Store the value of j and k along this slice in the common block
c           # comxyt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
c$$$        jcom = j
c$$$        kcom = k
c
c           # compute modifications fadd, gadd and hadd to fluxes along
c           # this slice:
c
            call flux3(1,maxm,meqn,maux,mbc,mx,
     &                 q1d,dtdx1d,dtdy,dtdz,aux1,aux2,aux3,
     &                 faddm,faddp,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),

     &                 work(i0bmcqxxm),work(i0bmcqxxp),
     &                 work(i0bpcqxxm),work(i0bpcqxxp),
     &                 work(i0cmcqxxm),work(i0cmcqxxp),
     &                 work(i0cpcqxxm),work(i0cpcqxxp),

     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cflgrid = dmax1(cflgrid,cfl1d)
c
c        # update fluxes for use in AMR:
c
            do 25 i=1,mx+1
               do 25 m=1,meqn
               fm(m,i,j,k) = fm(m,i,j,k) + faddm(m,i)
               fp(m,i,j,k) = fp(m,i,j,k) + faddp(m,i)
c
               gm(m,i,j  ,k-1) = gm(m,i,j  ,k-1) + gadd(m,i,1,-1)
               gp(m,i,j  ,k-1) = gp(m,i,j  ,k-1) + gadd(m,i,1,-1)
               gm(m,i,j  ,k  ) = gm(m,i,j  ,k  ) + gadd(m,i,1, 0)
               gp(m,i,j  ,k  ) = gp(m,i,j  ,k  ) + gadd(m,i,1, 0)
               gm(m,i,j  ,k+1) = gm(m,i,j  ,k+1) + gadd(m,i,1, 1)
               gp(m,i,j  ,k+1) = gp(m,i,j  ,k+1) + gadd(m,i,1, 1)
c
               gm(m,i,j+1,k-1) = gm(m,i,j+1,k-1) + gadd(m,i,2,-1)
               gp(m,i,j+1,k-1) = gp(m,i,j+1,k-1) + gadd(m,i,2,-1)
               gm(m,i,j+1,k  ) = gm(m,i,j+1,k  ) + gadd(m,i,2, 0)
               gp(m,i,j+1,k  ) = gp(m,i,j+1,k  ) + gadd(m,i,2, 0)
               gm(m,i,j+1,k+1) = gm(m,i,j+1,k+1) + gadd(m,i,2, 1)
               gp(m,i,j+1,k+1) = gp(m,i,j+1,k+1) + gadd(m,i,2, 1)
c
               hm(m,i,j-1,k) = hm(m,i,j-1,k) + hadd(m,i,1,-1)
               hp(m,i,j-1,k) = hp(m,i,j-1,k) + hadd(m,i,1,-1)
               hm(m,i,j  ,k) = hm(m,i,j  ,k) + hadd(m,i,1, 0)
               hp(m,i,j  ,k) = hp(m,i,j  ,k) + hadd(m,i,1, 0)
               hm(m,i,j+1,k) = hm(m,i,j+1,k) + hadd(m,i,1, 1)
               hp(m,i,j+1,k) = hp(m,i,j+1,k) + hadd(m,i,1, 1)
c
               hm(m,i,j-1,k+1) = hm(m,i,j-1,k+1) + hadd(m,i,2,-1)
               hp(m,i,j-1,k+1) = hp(m,i,j-1,k+1) + hadd(m,i,2,-1)
               hm(m,i,j  ,k+1) = hm(m,i,j  ,k+1) + hadd(m,i,2, 0)
               hp(m,i,j  ,k+1) = hp(m,i,j  ,k+1) + hadd(m,i,2, 0)
               hm(m,i,j+1,k+1) = hm(m,i,j+1,k+1) + hadd(m,i,2, 1)
               hp(m,i,j+1,k+1) = hp(m,i,j+1,k+1) + hadd(m,i,2, 1)

   25          continue
   50    continue
c
c
c     # perform y sweeps
c     ==================
c
c
      do 100 k = 0, mz+1
         do 100 i = 0, mx+1
c
            do 70 j = 1-mbc, my+mbc
               do 70 m=1,meqn
c                 # copy data along a slice into 1d array:
                  q1d(m,j) = qold(m,i,j,k)
   70          continue
c
         if (mcapa.gt.0)  then
           do 71 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(mcapa,i,j,k)
   71      continue
         endif
c
         if (maux .gt. 0)  then
            do 72 j = 1-mbc, my+mbc
               do 72 ma=1,maux
                 aux1(ma,j,1) = aux(ma,i-1,j,k-1)
                 aux1(ma,j,2) = aux(ma,i,j,k-1)
                 aux1(ma,j,3) = aux(ma,i+1,j,k-1)
                 aux2(ma,j,1) = aux(ma,i-1,j,k)
                 aux2(ma,j,2) = aux(ma,i,j,k)
                 aux2(ma,j,3) = aux(ma,i+1,j,k)
                 aux3(ma,j,1) = aux(ma,i-1,j,k+1)
                 aux3(ma,j,2) = aux(ma,i,j,k+1)
                 aux3(ma,j,3) = aux(ma,i+1,j,k+1)
   72          continue
         endif
c
c           # Store the value of i and k along this slice in the common block
c           # comxyzt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
c$$$            icom = i
c$$$            kcom = k
c
c           # compute modifications fadd, gadd and hadd to fluxes along this
c           # slice:
c
            call flux3(2,maxm,meqn,maux,mbc,my,
     &                 q1d,dtdy1d,dtdz,dtdx,aux1,aux2,aux3,
     &                 faddm,faddp,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),

     &                 work(i0bmcqxxm),work(i0bmcqxxp),
     &                 work(i0bpcqxxm),work(i0bpcqxxp),
     &                 work(i0cmcqxxm),work(i0cmcqxxp),
     &                 work(i0cpcqxxm),work(i0cpcqxxp),

     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cflgrid = dmax1(cflgrid,cfl1d)
c
c        # update fluxes for use in AMR:
c
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the g-fluxes
c           # gadd - modifies the h-fluxes
c           # hadd - modifies the f-fluxes
c
            do 75 j=1,my+1
               do 75 m=1,meqn
               gm(m,i,j,k) = gm(m,i,j,k) + faddm(m,j)
               gp(m,i,j,k) = gp(m,i,j,k) + faddp(m,j)
c
               hm(m,i-1,j,k) = hm(m,i-1,j,k) + gadd(m,j,1,-1)
               hp(m,i-1,j,k) = hp(m,i-1,j,k) + gadd(m,j,1,-1)
               hm(m,i  ,j,k) = hm(m,i  ,j,k) + gadd(m,j,1, 0)
               hp(m,i  ,j,k) = hp(m,i  ,j,k) + gadd(m,j,1, 0)
               hm(m,i+1,j,k) = hm(m,i+1,j,k) + gadd(m,j,1, 1)
               hp(m,i+1,j,k) = hp(m,i+1,j,k) + gadd(m,j,1, 1)
c
               hm(m,i-1,j,k+1) = hm(m,i-1,j,k+1) + gadd(m,j,2,-1)
               hp(m,i-1,j,k+1) = hp(m,i-1,j,k+1) + gadd(m,j,2,-1)
               hm(m,i  ,j,k+1) = hm(m,i  ,j,k+1) + gadd(m,j,2, 0)
               hp(m,i  ,j,k+1) = hp(m,i  ,j,k+1) + gadd(m,j,2, 0)
               hm(m,i+1,j,k+1) = hm(m,i+1,j,k+1) + gadd(m,j,2, 1)
               hp(m,i+1,j,k+1) = hp(m,i+1,j,k+1) + gadd(m,j,2, 1)
c
               fm(m,i  ,j,k-1) = fm(m,i  ,j,k-1) + hadd(m,j,1,-1)
               fp(m,i  ,j,k-1) = fp(m,i  ,j,k-1) + hadd(m,j,1,-1)
               fm(m,i  ,j,k  ) = fm(m,i  ,j,k  ) + hadd(m,j,1, 0)
               fp(m,i  ,j,k  ) = fp(m,i  ,j,k  ) + hadd(m,j,1, 0)
               fm(m,i  ,j,k+1) = fm(m,i  ,j,k+1) + hadd(m,j,1, 1)
               fp(m,i  ,j,k+1) = fp(m,i  ,j,k+1) + hadd(m,j,1, 1)
c
               fm(m,i+1,j,k-1) = fm(m,i+1,j,k-1) + hadd(m,j,2,-1)
               fp(m,i+1,j,k-1) = fp(m,i+1,j,k-1) + hadd(m,j,2,-1)
               fm(m,i+1,j,k  ) = fm(m,i+1,j,k  ) + hadd(m,j,2, 0)
               fp(m,i+1,j,k  ) = fp(m,i+1,j,k  ) + hadd(m,j,2, 0)
               fm(m,i+1,j,k+1) = fm(m,i+1,j,k+1) + hadd(m,j,2, 1)
               fp(m,i+1,j,k+1) = fp(m,i+1,j,k+1) + hadd(m,j,2, 1)
c

   75          continue
c
  100    continue
c
c
c
c     # perform z sweeps
c     ==================
c
c
      do 150 j = 0, my+1
         do 150 i = 0, mx+1
c
            do 110 k = 1-mbc, mz+mbc
               do 110 m=1,meqn
c                 # copy data along a slice into 1d array:
                  q1d(m,k) = qold(m,i,j,k)
 110           continue
c
         if (mcapa.gt.0)  then
           do 130 k = 1-mbc, mz+mbc
               dtdz1d(k) = dtdz / aux(mcapa,i,j,k)
 130       continue
         endif
c
         if (maux .gt. 0)  then
            do 131 k = 1-mbc, mz+mbc
               do 131 ma=1,maux
                 aux1(ma,k,1) = aux(ma,i-1,j-1,k)
                 aux1(ma,k,2) = aux(ma,i-1,j,k)
                 aux1(ma,k,3) = aux(ma,i-1,j+1,k)
                 aux2(ma,k,1) = aux(ma,i,j-1,k)
                 aux2(ma,k,2) = aux(ma,i,j,k)
                 aux2(ma,k,3) = aux(ma,i,j+1,k)
                 aux3(ma,k,1) = aux(ma,i+1,j-1,k)
                 aux3(ma,k,2) = aux(ma,i+1,j,k)
                 aux3(ma,k,3) = aux(ma,i+1,j+1,k)
  131          continue
           endif
c
c           # Store the value of i and j along this slice in the common block
c           # comxyzt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
c$$$            icom = i
c$$$            jcom = j
c
c           # compute modifications fadd, gadd and hadd to fluxes along this
c           # slice:
c
            call flux3(3,maxm,meqn,maux,mbc,mz,
     &                 q1d,dtdz1d,dtdx,dtdy,aux1,aux2,aux3,
     &                 faddm,faddp,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),

     &                 work(i0bmcqxxm),work(i0bmcqxxp),
     &                 work(i0bpcqxxm),work(i0bpcqxxp),
     &                 work(i0cmcqxxm),work(i0cmcqxxp),
     &                 work(i0cpcqxxm),work(i0cpcqxxp),

     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cflgrid = dmax1(cflgrid,cfl1d)
c
c
c        # update fluxes for use in AMR:
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the h-fluxes
c           # gadd - modifies the f-fluxes
c           # hadd - modifies the g-fluxes
c
            do 125 k=1,mz+1
               do 125 m=1,meqn
               hm(m,i,j,k) = hm(m,i,j,k) + faddm(m,k)
               hp(m,i,j,k) = hp(m,i,j,k) + faddp(m,k)
c
               fm(m,i  ,j-1,k) = fm(m,i  ,j-1,k) + gadd(m,k,1,-1)
               fp(m,i  ,j-1,k) = fp(m,i  ,j-1,k) + gadd(m,k,1,-1)
               fm(m,i  ,j  ,k) = fm(m,i  ,j  ,k) + gadd(m,k,1, 0)
               fp(m,i  ,j  ,k) = fp(m,i  ,j  ,k) + gadd(m,k,1, 0)
               fm(m,i  ,j+1,k) = fm(m,i  ,j+1,k) + gadd(m,k,1, 1)
               fp(m,i  ,j+1,k) = fp(m,i  ,j+1,k) + gadd(m,k,1, 1)
c
               fm(m,i+1,j-1,k) = fm(m,i+1,j-1,k) + gadd(m,k,2,-1)
               fp(m,i+1,j-1,k) = fp(m,i+1,j-1,k) + gadd(m,k,2,-1)
               fm(m,i+1,j  ,k) = fm(m,i+1,j  ,k) + gadd(m,k,2, 0)
               fp(m,i+1,j  ,k) = fp(m,i+1,j  ,k) + gadd(m,k,2, 0)
               fm(m,i+1,j+1,k) = fm(m,i+1,j+1,k) + gadd(m,k,2, 1)
               fp(m,i+1,j+1,k) = fp(m,i+1,j+1,k) + gadd(m,k,2, 1)
c
               gm(m,i-1,j  ,k) = gm(m,i-1,j  ,k) + hadd(m,k,1,-1)
               gp(m,i-1,j  ,k) = gp(m,i-1,j  ,k) + hadd(m,k,1,-1)
               gm(m,i  ,j  ,k) = gm(m,i  ,j  ,k) + hadd(m,k,1, 0)
               gp(m,i  ,j  ,k) = gp(m,i  ,j  ,k) + hadd(m,k,1, 0)
               gm(m,i+1,j  ,k) = gm(m,i+1,j  ,k) + hadd(m,k,1, 1)
               gp(m,i+1,j  ,k) = gp(m,i+1,j  ,k) + hadd(m,k,1, 1)
c
               gm(m,i-1,j+1,k) = gm(m,i-1,j+1,k) + hadd(m,k,2,-1)
               gp(m,i-1,j+1,k) = gp(m,i-1,j+1,k) + hadd(m,k,2,-1)
               gm(m,i  ,j+1,k) = gm(m,i  ,j+1,k) + hadd(m,k,2, 0)
               gp(m,i  ,j+1,k) = gp(m,i  ,j+1,k) + hadd(m,k,2, 0)
               gm(m,i+1,j+1,k) = gm(m,i+1,j+1,k) + hadd(m,k,2, 1)
               gp(m,i+1,j+1,k) = gp(m,i+1,j+1,k) + hadd(m,k,2, 1)
c
  125          continue
  150    continue
c
c
      return
      end
