c
c -------------------------------------------------------------
c
      subroutine stepgrid_dimSplit(q,fm,fp,gm,gp,hm,hp,mitot,mjtot,
     &                  mktot,mbc,dt,dtnew,dx,dy,dz,
     &                  nvar,xlow,ylow,zlow,time,mptr,maux,aux)
c
c
c ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
c                 dimensionally split version of stepgrid
c    take a step in x, then y, then z. for now only godunov splitting
c    not strang splitting. 
c
c take a time step on a single grid. overwrite solution array q.
c A modified version of the clawpack routine step3 is used.
c
c return fluxes in fm,fp and gm,gp and hm,hp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot by mktot).

c even though dimensionally split, need to save all 3 fluxes for
c flux conservation on the outside.

c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow,ylow,zlow = lower left corner of enlarged grid (including ghost cells).
c dt         = incoming time step
c dx,dy,dz   = mesh widths for this grid
c dtnew      = return suggested new time step for this grids soln.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use amr_module
      implicit double precision (a-h,o-z)

      external rpn3

      dimension q(nvar,mitot,mjtot,mktot)
      dimension fm(nvar,mitot,mjtot,mktot),fp(nvar,mitot,mjtot,mktot)
      dimension gm(nvar,mitot,mjtot,mktot),gp(nvar,mitot,mjtot,mktot)
      dimension hm(nvar,mitot,mjtot,mktot),hp(nvar,mitot,mjtot,mktot)
      dimension aux(maux,mitot,mjtot,mktot)

      logical    debug,  dump
      data       debug/.false./,  dump/.false./

c
c
!--        if (dump .and. mptr .ne. 1) 
!--     1       call prettyprint(q,nvar,mitot,mjtot,mktot,outunit)
      if (dump) then
         write(outunit,*)" grid ", mptr
         do k = 1, mktot
         do j = 1, mjtot
         do i = 1, mitot
           write(outunit,545) i,j,k,(q(ivar,i,j,k),ivar=1,nvar)
 545       format(3i3,3x,5e30.20)
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


c
c 2/28/02 : Added call to b4step3.
      call b4step3(mbc,mx,my,mz,nvar,q,
     &             xlowmbc,ylowmbc,zlowmbc,dx,dy,dz,time,dt,maux,aux)

c
c
c
c
c    :::  # take one step on the conservation law: Godunov splitting :::
c
      call step3x(mbig,nvar,maux,mbc,mx,my,mz,
     &           q,aux,dx,dt,cflgrid,
     &           fm,fp,rpn3)
c
!$OMP CRITICAL (setcfl)
      cfl_level = dmax1(cfl_level,cflgrid)
!$OMP END CRITICAL (setcfl)
c
c       # update q with x fluxes first.
      dtdx = dt/dx
        if (mcapa.eq.0) then
         do 50 k=1,mktot
         do 50 j=1,mjtot
         do 50 i=mbc+1,mitot-mbc
         do 50 m=1,nvar
c
c            # no capa array.  Standard flux differencing:

            q(m,i,j,k) = q(m,i,j,k) - dtdx*(fm(m,i+1,j,k)-fp(m,i,j,k))
     
 50       continue
       else
         do 51 k=1,mktot
         do 51 j=1,mjtot
         do 51 i=mbc+1,mitot-mbc
         do 51 m=1,nvar

c            # with capa array.

          q(m,i,j,k) = q(m,i,j,k)
     &              - dtdx * (fm(m,i+1,j,k) - fp(m,i,j,k))
     &              / aux(mcapa,i,j,k)
 51       continue
       endif
c
c  second step: update in y
      call step3y(mbig,nvar,maux,mbc,mx,my,mz,
     &           q,aux,dy,dt,cflgrid,
     &           gm,gp,rpn3)
c
!$OMP CRITICAL (setcfl)
      cfl_level = dmax1(cfl_level,cflgrid)
!$OMP END CRITICAL (setcfl)
c
c       # update q with y fluxes next
      dtdy = dt/dy
        if (mcapa.eq.0) then
         do 52 k=1,mktot
         do 52 j=mbc+1,mjtot-mbc
         do 52 i=mbc+1,mitot-mbc
         do 52 m=1,nvar
c
c            # no capa array.  Standard flux differencing:

            q(m,i,j,k) = q(m,i,j,k) - dtdy*(gm(m,i,j+1,k)-gp(m,i,j,k))
     
 52       continue
       else
         do 53 k=1,mktot
         do 53 j=mbc+1,mjtot-mbc
         do 53 i=mbc+1,mitot-mbc
         do 53 m=1,nvar

c            # with capa array.

          q(m,i,j,k) = q(m,i,j,k)
     &              - dtdy * (gm(m,i,j+1,k) - gp(m,i,j,k))
     &              / aux(mcapa,i,j,k)
 53       continue
       endif
c
c  third step: update in z
      call step3z(mbig,nvar,maux,mbc,mx,my,mz,
     &           q,aux,dz,dt,cflgrid,
     &           hm,hp,rpn3)
c
!$OMP CRITICAL (setcfl)
      cfl_level = dmax1(cfl_level,cflgrid)
!$OMP END CRITICAL (setcfl)
c
c       # update q with z fluxes last
      dtdz = dt/dz
        if (mcapa.eq.0) then
         do 54 k=mbc+1,mktot-mbc
         do 54 j=mbc+1,mjtot-mbc
         do 54 i=mbc+1,mitot-mbc
         do 54 m=1,nvar
c
c            # no capa array.  Standard flux differencing:

            q(m,i,j,k) = q(m,i,j,k) - dtdz*(hm(m,i,j,k+1)-hp(m,i,j,k))
     
 54       continue
       else
         do 55 k=mbc+1,mktot-mbc
         do 55 j=mbc+1,mjtot-mbc
         do 55 i=mbc+1,mitot-mbc
         do 55 m=1,nvar

c            # with capa array.

          q(m,i,j,k) = q(m,i,j,k)
     &              - dtdz * (hm(m,i,j,k+1) - hp(m,i,j,k))
     &              / aux(mcapa,i,j,k)
 55       continue
       endif
c
c     ::: end of Godunov dimensional splitting :::
c
      if (method(5).eq.1) then
c        # with source term:   use Godunov splitting
c         call src3(q,mitot,mjtot,mktot,nvar,aux,maux,time,dt)
         call src3(nvar,mbc,mx,my,mz,xlowmbc,ylowmbc,
     &             zlowmbc,dx,dy,dz,q,maux,aux,time,dt)
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
               write(dbugunit,831) i,j,k,fm(1,i,j,k),fp(1,i,j,k),
     .                           gm(1,i,j,k),gp(1,i,j,k),
     .                           hm(1,i,j,k),hp(1,i,j,k)
             do 830 m = 2, meqn
              write(dbugunit,832) fm(m,i,j,k),fp(m,i,j,k),
     .              gm(m,i,j,k),gp(m,i,j,k),
     .              hm(m,i,j,k),hp(m,i,j,k)
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
        if (cfl_level .gt. 0.d0) then
            dtnew = dt*cfl/cfl_level
          else
c            # velocities are all zero on this grid so there's no
c            # time step restriction coming from this grid.
            dtnew = rinfinity
          endif
c

      if (cfl_level .gt. cflv1) then
            write(*,810) cfl_level
            write(outunit,810) cfl_level
  810       format('*** WARNING *** Courant number  =', d12.4,
     &              '  is larger than cflv(1) ')
            endif

!--        if (dump .and. mptr .ne. 1) 
!--     1     call prettyprint(q,nvar,mitot,mjtot,mktot,outunit)
      if (dump) then
         write(outunit,*)" after time step on grid ", mptr
         do k = 1, mktot
         do j = 1, mjtot
         do i = 1, mitot
           write(outunit,545) i,j,k,(q(ivar,i,j,k),ivar=1,nvar)
         end do
         end do
         end do
      endif

      return
      end

