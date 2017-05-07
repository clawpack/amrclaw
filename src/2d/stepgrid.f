c
c -------------------------------------------------------------
c
!> Take a time step on a single grid **mptr** and overwrite solution array **q**. 
!! A modified version of the clawpack routine step2 is used.
!! 
!! Return new solution **q** as well as fluxes in fm,fp and gm,gp.
!! Patch has room for ghost cells (mbc of them) around the grid.
!! Everything is the enlarged size (**mitot** by **mjtot**).
!! 
!! \param[in] mbc number of ghost cells  (= lwidth)
!! \param[in] mptr grid number (for debugging)
!! \param[in] xlow left edge of enlarged grid (including ghost cells).
!! \param[in] ylow lower edge of enlarged grid (including ghost cells).
!! \param[in] dt incoming time step
!! \param[in] dx mesh size in x direction for this grid
!! \param[in] dx mesh size in y direction for this grid
!! \param[in,out] q solution array
!! \param[out] dtnew  return suggested new time step for this grid's soln.
!! \param[out] fm fluxes on the left side of each vertical edge
!! \param[out] fp fluxes on the right side of each vertical edge
!! \param[out] gm fluxes on the lower side of each horizontal edge
!! \param[out] gp fluxes on the upper side of each horizontal edge
      subroutine stepgrid(q,fm,fp,gm,gp,mitot,mjtot,mbc,dt,dtnew,dx,dy,
     &                  nvar,xlow,ylow,time,mptr,maux,aux)
c
c          
c ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
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
      external rpn2,rpt2

      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      parameter (msize=max1d+4)
      parameter (mwork=msize*(maxvar*maxvar + 13*maxvar + 3*maxaux +2))

      dimension q(nvar,mitot,mjtot)
      dimension fp(nvar,mitot,mjtot),gp(nvar,mitot,mjtot)
      dimension fm(nvar,mitot,mjtot),gm(nvar,mitot,mjtot)
      dimension aux(maux,mitot,mjtot)
c 	  dimension work(mwork)

      logical    debug,  dump
      data       debug/.false./,  dump/.false./

c
c     # set tcom = time.  This is in the common block comxyt that could
c     # be included in the Riemann solver, for example, if t is explicitly
c     # needed there.

      tcom = time

      if (dump) then
         write(outunit,*) "dumping grid ",mptr," at time ",time
         do i = 1, mitot
         do j = 1, mjtot
            write(outunit,545) i,j,(q(ivar,i,j),ivar=1,nvar) 
c    .                  ,(aux(ivar,i,j),ivar=1,maux)
 545        format(2i4,5e15.7)
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

c     # method(2:7) and mthlim
c     #    are set in the amr2ez file (read by amr)
c
      method(1) = 0
c
c
c     # fluxes initialized in step2
c
C       mwork0 = (maxm+2*mbc)*(12*meqn + mwaves + meqn*mwaves + 2) 
c
C       if (mwork .lt. mwork0) then
C          write(outunit,*) 'CLAW2 ERROR... mwork must be increased to ',
C      &               mwork0
C          write(*      ,*) 'CLAW2 ERROR... mwork must be increased to ',
C      &               mwork0
C          stop
C       endif
c  
c
c     # partition work array into pieces needed for local storage in 
c     # step2 routine. Find starting index of each piece:
c
C       i0faddm = 1
C       i0faddp = i0faddm + (maxm+2*mbc)*meqn
C       i0gaddm = i0faddp + (maxm+2*mbc)*meqn
C       i0gaddp = i0gaddm + 2*(maxm+2*mbc)*meqn
C       i0q1d   = i0gaddp + 2*(maxm+2*mbc)*meqn 
C       i0dtdx1 = i0q1d + (maxm+2*mbc)*meqn  
C       i0dtdy1 = i0dtdx1 + (maxm+2*mbc)
C       i0aux1 = i0dtdy1 + (maxm+2*mbc)
C       i0aux2 = i0aux1 + (maxm+2*mbc)*maux
C       i0aux3 = i0aux2 + (maxm+2*mbc)*maux
c
c
C       i0next = i0aux3 + (maxm+2*mbc)*maux    !# next free space
C       mused  = i0next - 1                    !# space already used
C       mwork1 = mwork - mused              !# remaining space (passed to step2)

c
c
      call b4step2(mbc,mx,my,nvar,q,
     &             xlowmbc,ylowmbc,dx,dy,time,dt,maux,aux)
c
c
c     # take one step on the conservation law:
c
      call step2(mbig,nvar,maux,
     &           mbc,mx,my,
     &              q,aux,dx,dy,dt,cflgrid,
     &              fm,fp,gm,gp,rpn2,rpt2)
c
c
c       write(outunit,1001) mptr, node(nestlevel,mptr),cflgrid
c1001   format(' Courant # of grid', i4,
c    &        ' on level', i3, ' is  ', e10.3)
c

!$OMP  CRITICAL (cflm)

        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)

c
c       # update q
        dtdx = dt/dx
        dtdy = dt/dy
        do 50 j=mbc+1,mjtot-mbc
        do 50 i=mbc+1,mitot-mbc
        do 50 m=1,nvar
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:

           q(m,i,j) = q(m,i,j) 
     &           - dtdx * (fm(m,i+1,j) - fp(m,i,j)) 
     &           - dtdy * (gm(m,i,j+1) - gp(m,i,j)) 
         else
c            # with capa array.
           q(m,i,j) = q(m,i,j) 
     &          - (dtdx * (fm(m,i+1,j) - fp(m,i,j))
     &          +  dtdy * (gm(m,i,j+1) - gp(m,i,j))) / aux(mcapa,i,j)
         endif

 50      continue
c
c
      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting
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
c
       if (cflgrid .gt. 0.d0) then
           dtnew = dt*cfl/cflgrid
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


