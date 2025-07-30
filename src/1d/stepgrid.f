c
c -------------------------------------------------------------
c
      subroutine stepgrid(q,fm,fp,mitot,mbc,dt,dtnew,dx,
     &                  nvar,xlow,time,mptr,maux,aux)
c
c          
c ::::::::::::::::::: STEPGRID ::::::::::::::::::::::::::::::::::::
c take a time step on a single grid. overwrite solution array q. 
c A modified version of the clawpack routine step is used.
c
c return fluxes in fm and fp.
c patch has room for ghost cells (mbc of them) around the grid.
c everything is the enlarged size (mitot by mjtot).
c
c mbc       = number of ghost cells  (= lwidth)
c mptr      = grid number  (for debugging)
c xlow      =  left corner of enlarged grid (including ghost cells).
c dt        = incoming time step
c dx        = mesh width for this grid
c dtnew     = return suggested new time step for this grid's soln.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      use amr_module
      implicit double precision (a-h,o-z)
      external rp1

      common /comxyt/ dtcom,dxcom,tcom,icom

      dimension q(nvar,mitot)
      dimension fp(nvar,mitot)
      dimension fm(nvar,mitot)
      dimension aux(maux,mitot)

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
            write(outunit,545) i,(q(ivar,i),ivar=1,nvar)
c    .                  ,(aux(ivar,i),ivar=1,maux)
 545        format(2i4,5e15.7)
         end do
      endif
c
      meqn   = nvar
      mx = mitot - 2*mbc
      mbig = mx       !# size for 1d scratch array
      xlowmbc = xlow + mbc*dx

c     # method(2:7) and mthlim
c     #    are set in the amr2ez file (read by amr)
c
      method(1) = 0
c
c     This call has been moved out to advanc
c      call b4step1(mbc,mx,nvar,q,
c     &             xlowmbc,dx,time,dt,maux,aux)
c
c
c     # take one step on the conservation law:
c
      call step1(mbig,nvar,maux,
     &           mbc,mx,
     &              q,aux,dx,dt,cflgrid,
     &              fm,fp,rp1)
c
c
        write(outunit,1001) mptr, node(nestlevel,mptr),cflgrid
 1001   format(' Courant # of grid', i4,
     &        ' on level', i3, ' is  ', e10.3)
c

!$OMP  CRITICAL (cflm)

        cfl_level = dmax1(cfl_level,cflgrid)

!$OMP END CRITICAL (cflm)

c
c       # update q
        dtdx = dt/dx
        do 50 i=mbc+1,mitot-mbc
        do 50 m=1,nvar
         if (mcapa.eq.0) then
c
c            # no capa array.  Standard flux differencing:

           q(m,i) = q(m,i)
     &           - dtdx * (fm(m,i+1) - fp(m,i))
         else
c            # with capa array.
           q(m,i) = q(m,i)
     &          - (dtdx * (fm(m,i+1) - fp(m,i))) / aux(mcapa,i)
         endif

 50      continue
c
c
      if (method(4).eq.1) then
c        # with source term:   use Godunov splitting
         call src1(nvar,mbc,mx,xlowmbc,dx,
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
               write(dbugunit,831) i,fm(1,i),fp(1,i)
               do 830 m = 2, meqn
                  write(dbugunit,832) fm(m,i),fp(m,i)
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
            write(outunit,545) i,(q(ivar,i),ivar=1,nvar)
         end do
      endif
      return
      end


