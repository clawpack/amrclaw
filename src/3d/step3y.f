c
c     ==================================================================
      subroutine step3y(maxm,meqn,maux,mbc,mx,my,mz,
     &                 qold,aux,dy,dt,cflgrid,
     &                 gm,gp,rpn3)
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

      external rpn3

      dimension
     & qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)
      dimension q1d(meqn, 1-mbc:maxm+mbc)
c
      dimension gm(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      dimension gp(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
c
      dimension faddm(meqn, 1-mbc:maxm+mbc)
      dimension faddp(meqn, 1-mbc:maxm+mbc)
c
      dimension
     & aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)
      dimension aux2(maux, 1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
c
c
      cflgrid = 0.d0
      dtdy = dt/dy
c
      do 10 i=1-mbc,mx+mbc
         do 10 j=1-mbc,my+mbc
            do 10 k=1-mbc,mz+mbc
               do 10 m=1,meqn
                  gm(m,i,j,k) = 0.d0
                  gp(m,i,j,k) = 0.d0
  10            continue

c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
            dtdy1d(i) = dtdy
 5       continue
      endif
c
c
c     # perform y-sweeps (for Godunov split method, called from stepgrid_dimSplit)
    ! loop indices assume x sweep goes first, then y
c     ==================
c
      do 100 k = 1-mbc, mz+mbc
         do 100 i = 1, mx
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
                 aux2(ma,j) = aux(ma,i,j,k) 
   72          continue
         endif
c

c           # compute modifications fadd to fluxes along this slice:
c
            call flux3_dimSplit(2,maxm,meqn,maux,mbc,my,
     &                 q1d,dtdy1d,aux2,
     &                 faddm,faddp,cfl1d,
     &                 rpn3)
c
            cflgrid = dmax1(cflgrid,cfl1d)
c
c        # update fluxes for use in AMR:
c
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the g-fluxes
c
            do 75 j=1,my+1
               do 75 m=1,meqn
               gm(m,i,j,k) = gm(m,i,j,k) + faddm(m,j)
               gp(m,i,j,k) = gp(m,i,j,k) + faddp(m,j)
   75          continue
c
  100    continue
c
      return
      end
