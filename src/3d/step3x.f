c
c     ==================================================================
      subroutine step3x(maxm,meqn,maux,mbc,mx,my,mz,
     &                 qold,aux,dx,dt,cflgrid,
     &                 fm,fp,rpn3)
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
      dimension fm(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      dimension fp(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
c
      dimension faddm(meqn, 1-mbc:maxm+mbc)
      dimension faddp(meqn, 1-mbc:maxm+mbc)
c
      dimension
     &      aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)
      dimension aux2(maux, 1-mbc:maxm+mbc) 

      dimension dtdx1d(1-mbc:maxm+mbc)
c
      cflgrid = 0.d0
      dtdx = dt/dx
c
      do 10 i=1-mbc,mx+mbc
         do 10 j=1-mbc,my+mbc
            do 10 k=1-mbc,mz+mbc
               do 10 m=1,meqn
                  fm(m,i,j,k) = 0.d0
                  fp(m,i,j,k) = 0.d0
  10            continue

c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
 5       continue
      endif
c
c
c     # perform x-sweeps (for Godunov split method, called from stepgrid_dimSplit)
    ! loop indices assume x sweep goes first
c     ==================
c
      do 50 k = 1-mbc,mz+mbc
         do 50 j = 1-mbc,my+mbc
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
                 aux2(ma,i) = aux(ma,i,j,k)
   22          continue
           endif
c
c
c           # compute modifications fadd to fluxes along this slice:
c
            call flux3_dimSplit(1,maxm,meqn,maux,mbc,mx,
     &                 q1d,dtdx1d,aux2,
     &                 faddm,faddp,cfl1d,
     &                 rpn3)
c
            cflgrid = dmax1(cflgrid,cfl1d)
c
c        # update fluxes for use in AMR:
c
            do 25 i=1,mx+1
               do 25 m=1,meqn
               fm(m,i,j,k) = fm(m,i,j,k) + faddm(m,i)
               fp(m,i,j,k) = fp(m,i,j,k) + faddp(m,i)
   25          continue
   50    continue
c
      return
      end
