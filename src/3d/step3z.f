c
c     ==================================================================
      subroutine step3z(maxm,meqn,maux,mbc,mx,my,mz,
     &                 qold,aux,dz,dt,cflgrid,
     &                 hm,hp,rpn3)
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
      dimension hm(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      dimension hp(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
c
      dimension faddm(meqn, 1-mbc:maxm+mbc)
      dimension faddp(meqn, 1-mbc:maxm+mbc)
c
      dimension
     & aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc)
      dimension aux2(maux, 1-mbc:maxm+mbc)
      dimension dtdz1d(1-mbc:maxm+mbc)
c
      cflgrid = 0.d0
      dtdz = dt/dz
c
      do 10 i=1-mbc,mx+mbc
         do 10 j=1-mbc,my+mbc
            do 10 k=1-mbc,mz+mbc
               do 10 m=1,meqn
                  hm(m,i,j,k) = 0.d0
                  hp(m,i,j,k) = 0.d0
 10            continue

c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
            dtdz1d(i) = dtdz
    5       continue
         endif
c
c
c     # perform z sweeps (for Godunov split method, called from stepgrid_dimSplit)
    ! loop indices assume x sweep goes first, then y. z goes last
c     ==================
c
c
      do 150 j = 1, my
         do 150 i = 1, mx
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
                 aux2(ma,k) = aux(ma,i,j,k) 
  131          continue
           endif
c
c           # compute modifications fadd to fluxes along this slice:
c
            call flux3_dimSplit(3,maxm,meqn,maux,mbc,mz,
     &                 q1d,dtdz1d,aux2,
     &                 faddm,faddp,cfl1d,
     &                 rpn3)
c
            cflgrid = dmax1(cflgrid,cfl1d)
c
c
c        # update fluxes for use in AMR:
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the h-fluxes
c
            do 125 k=1,mz+1
               do 125 m=1,meqn
               hm(m,i,j,k) = hm(m,i,j,k) + faddm(m,k)
               hp(m,i,j,k) = hp(m,i,j,k) + faddp(m,k)
  125          continue
  150    continue
c
c
      return
      end
