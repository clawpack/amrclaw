c
c -------------------------------------------------------------
c
      subroutine stepgrid_dimSplit(q,fm,fp,gm,gp,mitot,mjtot,mbc,
     &                             dt,dtnew,dx,dy,
     &                             nvar,xlow,ylow,time,mptr,maux,aux)
c
c          
c ::::::::::::::::::: STEPGRID_DIMSPLIT ::::::::::::::::::::::::::::::::::::
c                 dimensionally split version of stepgrid
c    take a step in x, then y. for now only godunov splitting
c    not strang splitting. 
c
c take a time step on a single grid. overwrite solution array q. 
c A modified version of the clawpack routine step2 is used.
c
c return fluxes in fm,fp and gm,gp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot).
c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow,ylow = lower left corner of enlarged grid (including ghost cells).
c dt         = incoming time step
c dx,dy      = mesh widths for this grid
c dtnew      = return suggested new time step for this grid's soln.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use amr_module
      implicit double precision (a-h,o-z)
      external rpn2

      parameter (msize=max1d+4)
      parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

      dimension q(nvar,mitot,mjtot)
      dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
      dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
C     dimension work(mwork)

      logical    debug,  dump
      data       debug/.false./,  dump/.false./

c
      if (dump) then
         write(outunit,*) "dumping grid ",mptr," at time ",time
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar),
     .                  (aux(ivar,i,j),ivar=1,maux)
 545        format(2i4,4e15.7)
         end do
         end do
      endif
c
      meqn   = nvar
      mx = mitot - 2*mbc
      my = mjtot - 2*mbc
      maxm = max(mx,my)       !# size for 1d scratch array
      mbig = maxm
      xlowmbc = xlow + mbc*dx
      ylowmbc = ylow + mbc*dy

c
c    old work array code now removed.
c
c
      call b4step2(mbc,mx,my,nvar,q,
     &             xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)
c
c
c     # take one step on the conservation law:  Godunov splitting
c
c                                   1.  First step in x direction:

      call step2x(mbig,nvar,maux, mbc,mx,my,q,aux,dx,dt,cflgrid,
     &              fm,fp,rpn2)
c
!$OMP  CRITICAL (cflm)

        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)

c       # update q with x fluxes first. all rows in y  get updated
        dtdx = dt/dx
        do 50 j=1,mjtot
        do 50 i=mbc+1,mitot-mbc
        do 50 m=1,nvar
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:

           q(m,i,j) = q(m,i,j)  - dtdx * (fm(m,i+1,j) - fp(m,i,j)) 

         else
c            # with capa array.
           q(m,i,j) = q(m,i,j) 
     &          - dtdx * (fm(m,i+1,j) - fp(m,i,j)) / aux(mcapa,i,j)
         endif

 50      continue

c                                 2.  Second step in x direction:

      call step2y(mbig,nvar,maux, mbc,mx,my,q,aux,dy,dt,cflgrid,
     &              gm,gp,rpn2)
c
!$OMP  CRITICAL (cflm)

        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)
c
c       # update q
        dtdy = dt/dy
        do 51 j=mbc+1,mjtot-mbc
        do 51 i=mbc+1,mitot-mbc
        do 51 m=1,nvar
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:

           q(m,i,j) = q(m,i,j) - dtdy * (gm(m,i,j+1) - gp(m,i,j)) 
         else
c            # with capa array.
           q(m,i,j) = q(m,i,j) 
     &          - dtdy * (gm(m,i,j+1) - gp(m,i,j)) / aux(mcapa,i,j)
         endif

 51      continue
c
!--        write(outunit,1001) mptr, node(nestlevel,mptr),cflgrid
!-- 1001   format(' Courant # of grid', i4,
!--     &        ' on level', i3, ' is  ', e10.3)

c
      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting for this too
         call src2(nvar,mbc,mx,my,xlowmbc,ylowmbc,dx,dy,
     &             q,maux,aux,time,dt)
         endif
c
c
c
c     # output fluxes for debugging purposes:
      if (debug) then
         write(dbugunit,*)" fluxes for grid ",mptr
c        do 830 j = mbc+1, mjtot-1
            do 830 i = mbc+1, mitot-1
         do 830 j = mbc+1, mjtot-1
               write(dbugunit,831) i,j,fm(1,i,j),fp(1,i,j),
     .                             gm(1,i,j),gp(1,i,j)
               do 830 m = 2, meqn
                  write(dbugunit,832) fm(m,i,j),fp(m,i,j),
     .            gm(m,i,j),gp(m,i,j)
  831          format(2i4,4d16.6)
  832          format(8x,4d16.6)
  830    continue
      endif

c
c
c For variable time stepping, use max speed seen on this grid to 
c choose the allowable new time step dtnew.  This will later be 
c compared to values seen on other grids.
c  NB: orig code tested cflgrid, but with dimensional splitting
c  there are 2 - one for each x and y steps, so changed 12/23/14 (mjb)
c  to use level
c
       if (cfl_level .gt. 0.d0) then
           dtnew = dt*cfl/cfl_level
         else
c          # velocities are all zero on this grid so there's no 
c          # time step restriction coming from this grid.
            dtnew = rinfinity
          endif

c     # give a warning if Courant number too large...
c
      if (cflgrid .gt. cflv1) then
            write(*,810) cflgrid
            write(outunit,810) cflgrid, cflv1
  810       format('*** WARNING *** Courant number  =', d12.4,
     &              '  is larger than input cfl_max = ', d12.4)
            endif
c
      if (dump) then
         write(outunit,*) "dumping grid ",mptr," after stepgrid"
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar)
         end do
         end do
      endif
      return
      end


