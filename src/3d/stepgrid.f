c
c -------------------------------------------------------------
c
      subroutine stepgrid(q,fm,fp,gm,gp,hm,hp,mitot,mjtot,mktot,
     &                  mbc,dt,dtnew,dx,dy,dz,
     &                  nvar,xlow,ylow,zlow,time,mptr,maux,aux)
c
c
c ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
c take a time step on a single grid. overwrite solution array q.
c A modified version of the clawpack routine step3 is used.
c
c return fluxes in fm,fp and gm,gp and hm,hp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot by mktot).
c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow,ylow,zlow = lower left corner of enlarged grid (including ghost cells).
c dt         = incoming time step
c dx,dy,dz   = mesh widths for this grid
c dtnew      = return suggested new time step for this grid's soln.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      implicit double precision (a-h,o-z)
      external rpn3,rpt3, rptt3

      include  "call.i"
      common/comxyzt/dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

      parameter (msize=max1d+4)

      parameter (mwork=msize*(46*maxvar + (maxvar+1)*maxwave
     &                 + 9*maxaux + 3))

      dimension q(mitot,mjtot,mktot,nvar)
      dimension fm(mitot,mjtot,mktot,nvar),fp(mitot,mjtot,mktot,nvar)
      dimension gm(mitot,mjtot,mktot,nvar),gp(mitot,mjtot,mktot,nvar)
      dimension hm(mitot,mjtot,mktot,nvar),hp(mitot,mjtot,mktot,nvar)
      dimension aux(mitot,mjtot,mktot,maux)
      dimension work(mwork)

      logical    debug,  dump
      data       debug/.false./,  dump/.false./
c
c     # set tcom = time.  This is in the common block comxyt that could
c     # be included in the Riemann solver, for example, if t is explicitly
c     # needed there.

      tcom = time

c
      if (dump) then
         write(outunit,*)" grid ", mptr
	 do k = nghost+1, mktot-nghost
	 do j = nghost+1, mjtot-nghost
	 do i = nghost+1, mitot-nghost
	    write(outunit,545) i,j,k,(q(i,j,k,ivar),ivar=1,nvar)
 545        format(3i3,3x,5e15.7)
	 end do
	 end do
	 end do
      endif
c
      meqn   = nvar
      mx = mitot - 2*mbc
      my = mjtot - 2*mbc
      mz = mktot - 2*mbc
      maxm = max(mx,my,mz)       !# size for 1d scratch array
      mbig = maxm
      xlowmbc = xlow + mbc*dx
      ylowmbc = ylow + mbc*dy
      zlowmbc = zlow + mbc*dz

      if (maxm+2*mbc > max1d + 4) then
         write(6,*) 'Increase size of max1d; ', maxm+2*mbc
         stop
      endif


c     # method(2), method(3), method(5), and mthlim
c     #    are set in the amr3ez.data file (read by amr)
c

      method(1) = 0
c
c
c     # partition work array into pieces needed for local storage in
c     # step3 routine. Find starting index of each piece:
c
      i0faddm = 1
      i0faddp = i0faddm   + (maxm+2*mbc)*meqn
      i0gadd  = i0faddp   + (maxm+2*mbc)*meqn
      i0hadd  = i0gadd  + 6*(maxm+2*mbc)*meqn
      i0q1d   = i0hadd  + 6*(maxm+2*mbc)*meqn
      i0dtdx1 = i0q1d     + (maxm+2*mbc)*meqn
      i0dtdy1 = i0dtdx1   + (maxm+2*mbc)
      i0dtdz1 = i0dtdy1   + (maxm+2*mbc)
      i0aux1 = i0dtdz1    + (maxm+2*mbc)
      i0aux2 = i0aux1     + (maxm+2*mbc)*maux*3
      i0aux3 = i0aux2     + (maxm+2*mbc)*maux*3
c
c
      i0next = i0aux3 + (maxm+2*mbc)*maux*3    !# next free space
      mused  = i0next - 1                    !# space already used
      mwork1 = mwork - mused              !# remaining space (passed to step3)

c 2/28/02 : Added call to b4step3.
      call b4step3(mx,my,mz,mbc,mx,my,mz,nvar,q,
     &             xlowmbc,ylowmbc,zlowmbc,dx,dy,dz,time,dt,maux,aux)

c
c
c
c
c     # take one step on the conservation law:
c
      call step3(mbig,mx,my,mz,nvar,maux,mbc,mx,my,mz,
     &		 q,aux,dx,dy,dz,dt,cflgrid,
     &           fm,fp,gm,gp,hm,hp,
     &           work(i0faddm),work(i0faddp),
     &           work(i0gadd),work(i0hadd),
     &           work(i0q1d),work(i0dtdx1),work(i0dtdy1),work(i0dtdz1),
     &		 work(i0aux1),work(i0aux2),work(i0aux3),
     &           work(i0next),mwork1,rpn3,rpt3, rptt3)
c
c
c       write(outunit,*) ' Courant # of grid ',mptr, '  is  ',cflgrid
c
	cflmax = dmax1(cflmax,cflgrid)
c
c       # update q
	dtdx = dt/dx
	dtdy = dt/dy
	dtdz = dt/dz
        if (mcapa.eq.0) then
	   do 50 m=1,nvar
	   do 50 i=mbc+1,mitot-mbc
	   do 50 j=mbc+1,mjtot-mbc
	   do 50 k=mbc+1,mktot-mbc
c
c            # no capa array.  Standard flux differencing:

	      q(i,j,k,m) = q(i,j,k,m)
     &	  	   - dtdx * (fm(i+1,j,k,m) - fp(i,j,k,m))
     &	  	   - dtdy * (gm(i,j+1,k,m) - gp(i,j,k,m))
     &	  	   - dtdz * (hm(i,j,k+1,m) - hp(i,j,k,m))
 50       continue
       else
	   do 51 m=1,nvar
	   do 51 i=mbc+1,mitot-mbc
	   do 51 j=mbc+1,mjtot-mbc
	   do 51 k=mbc+1,mktot-mbc

c            # with capa array.

	    q(i,j,k,m) = q(i,j,k,m)
     &		 - (dtdx * (fm(i+1,j,k,m) - fp(i,j,k,m))
     & 		 +  dtdy * (gm(i,j+1,k,m) - gp(i,j,k,m))
     &		 +  dtdz * (hm(i,j,k+1,m) - hp(i,j,k,m)))
     &         / aux(i,j,k,mcapa)
 51       continue
       endif

c
c
      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting
c         call src3(q,mitot,mjtot,mktot,nvar,aux,maux,time,dt)
         call src3(mx,my,mz,nvar,mbc,mx,my,mz,xlowmbc,
     &         ylowmbc,zlowmbc,dx,dy,dz,q,maux,aux,time,dt)
	 endif
c
c
c
c     # output fluxes for debugging purposes:
      if (debug) then
	 write(dbugunit,*)" fluxes for grid ",mptr
         do 830 i = mbc+1, mitot-1
            do 830 j = mbc+1, mjtot-1
            do 830 k = mbc+1, mktot-1
               write(dbugunit,831) i,j,k,fm(i,j,k,1),fp(i,j,k,1),
     .				   gm(i,j,k,1),gp(i,j,k,1),
     .				   hm(i,j,k,1),hp(i,j,k,1)
	       do 830 m = 2, meqn
		  write(dbugunit,832) fm(i,j,k,m),fp(i,j,k,m),
     .		  gm(i,j,k,m),gp(i,j,k,m),
     .		  hm(i,j,k,m),hp(i,j,k,m)
  831          format(3i4,6d16.6)
  832          format(8x,6d16.6)
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
c            # velocities are all zero on this grid so there's no
c            # time step restriction coming from this grid.
 	     dtnew = rinfinity
 	   endif
c

      if (cflgrid .gt. cflv1) then
            write(*,810) cflgrid
            write(outunit,810) cflgrid
  810       format('*** WARNING *** Courant number  =', d12.4,
     &              '  is larger than cflv(1) ')
            endif

      return
      end
