c
c ------------------------------------------------------------------
c
      subroutine bc3amr(val,aux,nrow,ncol,nfil,meqn,naux,
     1                  hx, hy, hz, level, time,
     2                  xleft,  xright,  yfront, yrear,
     3                  zbot, ztop,
     4                  xlower,ylower,zlower,
     5                  xupper,yupper,zupper,
     6                  xperiodic, yperiodic,zperiodic)


c
c
c :::::::::: BC3AMR ::::::::::::::::::::::::::::::::::::::::::::::;
c
c     Take a grid patch with mesh widths hx,hy,hz, of dimensions nrow by
c     ncol by nfil,  and set the values of any piece of
c     of the patch which extends outside the physical domain
c     using the boundary conditions.
c     ------------------------------------------------
c     # Standard boundary condition choices for amr3ez in clawpack
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (front), 4 (rear),
c                         5 (bottom) 6 (top)
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  or 4'th (for k = 5,6) component of q.
c     ------------------------------------------------
c
c     The corners of the grid patch are at
c        (xleft,yfront,zbot)  --  lower front left corner
c        (xright,yrear,ztop) --  upper rear right corner
c
c     The physical domain itself is a rectangular parallelopiped bounded by
c        (xlower,ylower,zlower)  -- lower front left corner
c        (xupper,yupper,zupper)  -- upper rear right corner
c
c     the picture is the following:
c
c                            __________________________(xupper,yupper,zupper)
c                           /                         /|
c                          /                         / |
c                         /                         /  |
c                        /_________________________/   |
c                        |                         |   |
c                        |                         |   |
c                     ___|_____(xright,yrear,ztop) |   |
c                    /___|____/|                   |   |
c                    |   |    ||                   |   |
c                    |   |    ||                   |   |
c                    |   |    ||                   |   |
c                    |___|____|/                   |   |
c  (xleft,yfront,zbot)   |                         |  /
c                        |                         | /
c                        |_________________________|/
c  (xlower,ylower,zlower)
c
c     Any cells that lie outside the physical domain are ghost cells whose
c     values should be set in this routine.  This is tested for by comparing
c     xleft with xlower to see if values need to be set at the left, as in
c     the figure above, and similarly at the other boundaries.
c
c     Patches are guaranteed to have at least 1 row of cells filled
c     with interior values so it is possible to  extrapolate.
c     Fix trimbd if you want more than 1 row pre-set.
c
c     Make sure the order the boundaries are specified is correct
c     so that diagonal corner cells are also properly taken care of.
c
c     Periodic boundaries are set before calling this routine, so if the
c     domain is periodic in one direction only you
c     can safely extrapolate in the other direction.
c
c     Don't overwrite ghost cells in periodic directions!
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

      use amr_module, only:  mthbc
      implicit none

      real*8  val(meqn,nrow,ncol,nfil), aux(naux,nrow,ncol,nfil)
      real*8  ubar, vbar, wbar
      logical xperiodic, yperiodic, zperiodic
      integer nrow,ncol,nfil,meqn,naux,level
      real*8  hx,hy,hz,time,qtrue
      real*8  hxmarg,hymarg,hzmarg
      real*8  xleft,xright,yfront,yrear,zbot,ztop
      real*8  xcell,ycell,zcell
      real*8  xlower,ylower,zlower
      real*8  xupper,yupper,zupper
      integer i,j,k,m,nxl,nxr,ibeg,nyf,nyr,jbeg,nzb,nzt,kbeg
      
      common /cparam/ ubar,vbar,wbar

      hxmarg = hx*.01d0
      hymarg = hy*.01d0
      hzmarg = hz*.01d0

      if (xperiodic .and. yperiodic .and. zperiodic) go to 699
c
c
c-------------------------------------------------------
c     # xlower boundary:
c-------------------------------------------------------
      if (xleft .ge. xlower-hxmarg) then
c        # not a physical boundary -- ghost cells lie within another
c        # grid and values are set elsewhere in amr code.
         go to 199
         endif
c
c     # number of ghost cells lying outside physical domain:
      nxl = (xlower+hxmarg-xleft)/hx     
c
      if (nxl > 2) then
          write(6,*) '*** unexpected value nxl = ',nxl
          stop
          endif
c
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
         do k = 1,nfil
            zcell = zbot + (k-0.5d0)*hz
            do j = 1,ncol
               ycell = yfront + (j-0.5d0)*hy
               do i=1,nxl 
                  xcell = xleft + (i-0.5d0)*hx
                  val(1,i,j,k) = qtrue(xcell,ycell,zcell,time)
               enddo
            enddo
         enddo
      go to 199
c
  110 continue
c     # zero-order extrapolation:
         do 115 k = 1,nfil
          do 115 j = 1,ncol
           do 115 i=1,nxl
            do 115 m=1,meqn
               val(m,i,j,k) = val(m,nxl+1,j,k)
  115       continue
      go to 199

  120 continue
c     # periodic:   handled elsewhere in amr
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 k = 1,nfil
       do 135 j = 1,ncol
        do 135 i=1,nxl
         do 135 m=1,meqn
               val(m,i,j,k) = val(m,2*nxl+1-i,j,k)
  135    continue
c     # negate the normal velocity:
      do 136 i=1,nxl
         do 136 j = 1,ncol
            do 136 k = 1,nfil
               val(2,i,j,k) = -val(2,i,j,k)
  136       continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # xupper boundary:
c-------------------------------------------------------
      if (xright .le. xupper+hxmarg) then
c        # not a physical boundary -- ghost cells lie within another
c        # grid and values are set elsewhere in amr code.
         go to 299
         endif
c
c     # number of ghost cells lying outside physical domain:
      nxr = (xright - xupper + hxmarg)/hx
      ibeg = max0(nrow-nxr+1, 1)
c
      go to (200,210,220,230) mthbc(2)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*)
     &   '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2amr'
      stop
      go to 299

  210 continue
c     # linear-order extrapolation:
       do 215 k = 1,nfil
        do 215 j = 1,ncol
         do 215 i=ibeg,nrow
          do 215 m=1,meqn
            val(m,i,j,k) = 2.d0*val(m,ibeg-1,j,k) - val(m,ibeg-2,j,k)
  215     continue
      go to 299

  220 continue
c     # periodic:   handled elsewhere in amr
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 k = 1,nfil
       do 235 j = 1,ncol
        do 235 i=ibeg,nrow
         do 235 m=1,meqn
            val(m,i,j,k) = val(m,2*ibeg-1-i,j,k)
  235    continue
c     # negate the normal velocity:
       do 236 k = 1,nfil
       do 236 j = 1,ncol
       do 236 i=ibeg,nrow
          val(2,i,j,k) = -val(2,i,j,k)
  236  continue
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # ylower boundary:
c-------------------------------------------------------
      if (yfront .ge. ylower-hymarg) then
c        # not a physical boundary -- ghost cells lie within another
c        # grid and values are set elsewhere in amr code.
         go to 399
         endif
c
c     # number of ghost cells lying outside physical domain:
      nyf = (ylower+hymarg-yfront)/hy
c
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
         do k = 1,nfil
            zcell = zbot + (k-0.5d0)*hz
            do i = 1,nrow
               xcell = xleft + (i-0.5d0)*hx
               do j=1,nyf
                   ycell = yfront + (j-.5d0)*hy
                   val(1,i,j,k) = qtrue(xcell,ycell,zcell,time)
               enddo
            enddo
         enddo
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 k = 1,nfil
       do 315 j=1,nyf
        do 315 i=1,nrow
         do 315 m=1,meqn
            val(m,i,j,k) = val(m,i,nyf+1,k)
  315 continue
      go to 399

  320 continue
c     # periodic:   handled elsewhere in amr
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 k = 1,nfil
       do 335 j=1,nyf
        do 335 i=1,nrow
         do 335 m=1,meqn
            val(m,i,j,k) =  val(m,i,2*nyf+1-j,k)
  335    continue
c     # negate the normal velocity:
      do 336 k = 1,nfil
       do 336 j=1,nyf
        do 336 i=1,nrow
           val(3,i,j,k) = -val(3,i,j,k)
  336   continue
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # yupper boundary:
c-------------------------------------------------------
      if (yrear .le. yupper+hymarg) then
c        # not a physical boundary -- ghost cells lie within another
c        # grid and values are set elsewhere in amr code.
         go to 499
         endif
c
c     # number of ghost cells lying outside physical domain:
      nyr = (yrear - yupper + hymarg)/hy
      jbeg = max0(ncol-nyr+1, 1)
c
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*)
     &   '*** ERROR *** mthbc(4)=0 and no BCs specified in bc3amr'
      stop
      go to 499

  410 continue
c     # linear  extrapolation:
       do 415 k = 1,nfil
        do 415 j=jbeg,ncol
         do 415 i=1,nrow
          do 415 m=1,meqn
            val(m,i,j,k) =  2.d0*val(m,i,jbeg-1,k) - val(m,i,jbeg-2,k)
  415     continue
      go to 499

  420 continue
c     # periodic:   handled elsewhere in amr
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 m=1,meqn
         do 435 j=jbeg,ncol
            do 435 i=1,nrow
               do 435 k = 1,nfil
                  val(m,i,j,k) =  val(m,i,2*jbeg-1-j,k)
  435          continue
c     # negate the normal velocity:
      do 436 j=jbeg,ncol
         do 436 i=1,nrow
            do 436 k = 1,nfil
               val(3,i,j,k) = -val(3,i,j,k)
  436       continue
      go to 499

  499 continue

c
c-------------------------------------------------------
c     # zlower boundary:
c-------------------------------------------------------
      if (zbot .ge. zlower-hzmarg) then
c        # not a physical boundary -- ghost cells lie within another
c        # grid and values are set elsewhere in amr code.
         go to 599
         endif
c
c     # number of ghost cells lying outside physical domain:
      nzb = (zlower+hzmarg-zbot)/hz
c
      go to (500,510,520,530) mthbc(5)+1
c
  500 continue
         do i = 1,nrow
            xcell = xleft + (i-0.5d0)*hx
            do j = 1,ncol
               ycell = yfront + (j-0.5d0)*hy
               do k=1,nzb
                  zcell = zbot + (k-0.5d0)*hz
                  val(1,i,j,k) = qtrue(xcell,ycell,zcell,time)
               enddo
            enddo
         enddo
      go to 599
c
  510 continue
c     # zero-order extrapolation:
      do 515 k=1,nzb
       do 515 j = 1,ncol
        do 515 i=1,nrow
         do 515 m=1,meqn
            val(m,i,j,k) = val(m,i,j,nzb+1)
  515    continue
      go to 599

  520 continue
c     # periodic:   handled elsewhere in amr
      go to 599

  530 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 535 k=1,nzb
       do 535 j = 1,ncol
        do 535 i=1,nrow
         do 535 m=1,meqn
            val(m,i,j,k) =  val(m,i,j,2*nzb+1-k)
  535    continue
c     # negate the normal velocity:
      do 536 k = 1,nzb
       do 536 j = 1,ncol
        do 536 i = 1,nrow
           val(4,i,j,k) = -val(4,i,j,k)
  536   continue
      go to 599

  599 continue
c
c-------------------------------------------------------
c     # zupper boundary:
c-------------------------------------------------------
      if (ztop .le. zupper+hzmarg) then
c        # not a physical boundary -- ghost cells lie within another
c        # grid and values are set elsewhere in amr code.
         go to 699
         endif
c
c     # number of ghost cells lying outside physical domain:
      nzt = (ztop - zupper + hzmarg)/hz
      kbeg = max0(nfil-nzt+1, 1)
c
      go to (600,610,620,630) mthbc(6)+1
c
  600 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*)
     &   '*** ERROR *** mthbc(6)=0 and no BCs specified in bc3amr'
      stop
      go to 699

  610 continue
c     # linear extrapolation:
      do 615 k = kbeg,nfil
       do 615 j = 1,ncol
        do 615 i = 1,nrow
         do 615 m = 1,meqn
           val(m,i,j,k) =  2.d0*val(m,i,j,kbeg-1) - val(m,i,j,kbeg-2)
  615    continue
      go to 699

  620 continue
c     # periodic:   handled elsewhere in amr
      go to 699

  630 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
       do 635 k = kbeg,nfil            
        do 635 j = 1,ncol
         do 635 i = 1,nrow
          do 635 m=1,meqn
             val(m,i,j,k) =  val(m,i,j,2*kbeg-1-k)
  635     continue
c     # negate the normal velocity:
      do 636 k=kbeg,nfil
       do 636 j = 1,ncol
        do 636 i=1,nrow
           val(4,i,j,k) = -val(4,i,j,k)
  636   continue
      go to 699

  699 continue

      return
      end
