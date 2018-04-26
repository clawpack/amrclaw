! :::::::::: bc2amr ::::::::::::::::::::::::::::::::::::::::::::::;
!> \callgraph
!! \callergraph
!!  Take a grid patch with mesh widths **hx**,**hy**, of dimensions **nrow** by
!!  **ncol**,  and set the values of any piece of
!!  of the patch which extends outside the physical domain 
!!  using the boundary conditions. 
!!
!!  ### Standard boundary condition choices for amr2ez in clawpack
!!
!!  At each boundary  k = 1 (left),  2 (right),  3 (bottom), 4 (top):
!!
!!  mthbc(k) =  
!!  * 0  for user-supplied BC's (must be inserted!)
!!  * 1  for zero-order extrapolation
!!  * 2  for periodic boundary conditions
!!  * 3  for solid walls, assuming this can be implemented
!!                   by reflecting the data about the boundary and then
!!                   negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
!!                   component of q.
!!  * 5  for sphere bcs (left half maps to right half of same side, and vice versa), as if domain folded in half
!!
!!  The corners of the grid patch are at 
!!     (xlo_patch,ylo_patch)  --  lower left corner
!!     (xhi_patch,yhi_patch) --  upper right corner
!!
!!  The physical domain itself is a rectangle bounded by
!!     (xlower,ylower)  -- lower left corner
!!     (xupper,yupper)  -- upper right corner
!!  
!   This figure below does not work with doxygen
!   the picture is the following: 
!  ____________________________________________________
! 
!                _____________________ (xupper,yupper)
!               |                     |  
!           ____|____ (xhi_patch,yhi_patch)   
!           |   |    |                |
!           |   |    |                |
!           |   |    |                |
!           |___|____|                |
!  (xlo_patch,ylo_patch) |            |
!               |                     |
!               |_____________________|   
!    (xlower,ylower)
!  ____________________________________________________
!!
!!
!>  Any cells that lie outside the physical domain are ghost cells whose
!!  values should be set in this routine.  This is tested for by comparing
!!  xlo_patch with xlower to see if values need to be set at the left
!   as in the figure above, 
!
!>  and similarly at the other boundaries.
!!  Patches are guaranteed to have at least 1 row of cells filled
!!  with interior values so it is possible to extrapolate. 
!!  Fix [trimbd()](@ref trimbd) if you want more than 1 row pre-set.
!!
!!  Make sure the order the boundaries are specified is correct
!!  so that diagonal corner cells are also properly taken care of.
!!
!!  Periodic boundaries are set before calling this routine, so if the
!!  domain is periodic in one direction only you
!!  can safely extrapolate in the other direction. 
!!
!!  Don't overwrite ghost cells in periodic directions!
!!
!! \param val data array for solution \f$q \f$ (cover the whole grid **msrc**)
!! \param aux data array for auxiliary variables 
!! \param nrow number of cells in *i* direction on this grid
!! \param ncol number of cells in *j* direction on this grid
!! \param meqn number of equations for the system
!! \param naux number of auxiliary variables
!! \param hx spacing (mesh size) in *i* direction
!! \param hy spacing (mesh size) in *j* direction
!! \param level AMR level of this grid
!! \param time setting ghost cell values at time **time**
!! \param xlo_patch left bound of the input grid
!! \param xhi_patch right bound of the input grid 
!! \param ylo_patch lower bound of the input grid 
!! \param yhi_patch upper bound of the input grid 
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;

c
c ------------------------------------------------------------------
c
      subroutine bc2amr(val,aux,nrow,ncol,meqn,naux,
     1                  hx, hy, level, time, 
     2                  xlo_patch, xhi_patch, ylo_patch,yhi_patch)
 
c
c

      use amr_module, only: mthbc,xlower,ylower,xupper,yupper
      use amr_module, only: xperdom,yperdom,spheredom

      implicit none

      real(CLAW_REAL)  val(meqn,nrow,ncol), aux(naux,nrow,ncol)
      integer nrow,ncol,meqn,naux,level
      real(CLAW_REAL)  hx,hy,time, hxmarg, hymarg
      real(CLAW_REAL)  xlo_patch,xhi_patch,ylo_patch,yhi_patch
      integer nxl,nxr,ibeg,nyb,nyt,jbeg,i,j,m
      

      hxmarg = hx*.01
      hymarg = hy*.01

      if (xperdom .and. (yperdom .or. spheredom)) go to 499
c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      if (xlo_patch .ge. xlower-hxmarg) then
c        # not a physical boundary -- no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 199
         endif
c
c     # number of grid cells from this patch lying outside physical domain:
      nxl = (xlower+hxmarg-xlo_patch)/hx
c
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) 
     &   '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2amr'
      stop
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do 115 j = 1,ncol
         do 115 i=1,nxl
            do 115 m=1,meqn
               val(m,i,j) = val(m,nxl+1,j)
  115       continue
      go to 199

  120 continue
c     # periodic:   handled elsewhere in amr
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 j = 1,ncol
         do 135 i=1,nxl
            do 135 m=1,meqn
               val(m,i,j) = val(m,2*nxl+1-i,j)
  135       continue
c     # negate the normal velocity:
      do 136 j = 1,ncol
         do 136 i=1,nxl
            val(2,i,j) = -val(2,i,j)
  136    continue
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      if (xhi_patch .le. xupper+hxmarg) then
c        # not a physical boundary --  no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 299
         endif
c
c     # number of grid cells lying outside physical domain:
      nxr = (xhi_patch - xupper + hxmarg)/hx
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
c     # zero-order extrapolation:
      do 215 j = 1,ncol
         do 215 i=ibeg,nrow
            do 215 m=1,meqn
               val(m,i,j) = val(m,ibeg-1,j)
  215       continue
      go to 299
  
  220 continue
c     # periodic:   handled elsewhere in amr
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 j = 1,ncol
         do 235 i=ibeg,nrow
            do 235 m=1,meqn
               val(m,i,j) = val(m,2*ibeg-1-i,j)
  235       continue
c     # negate the normal velocity:
      do 236 j = 1,ncol
         do 236 i=ibeg,nrow
            val(2,i,j) = -val(2,i,j)
  236    continue
      go to 299
   
  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      if (ylo_patch .ge. ylower-hymarg) then
c        # not a physical boundary -- no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 399
         endif
c
c     # number of grid cells lying outside physical domain:
      nyb = (ylower+hymarg-ylo_patch)/hy
c
      go to (300,310,320,330) mthbc(3)+1
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) 
     &   '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2amr'
      stop
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 j=1,nyb
         do 315 i=1,nrow
            do 315 m=1,meqn
                val(m,i,j) = val(m,i,nyb+1)
  315       continue
      go to 399

  320 continue
c     # periodic:   handled elsewhere in amr
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 j=1,nyb
         do 335 i=1,nrow
            do 335 m=1,meqn
               val(m,i,j) =  val(m,i,2*nyb+1-j)
  335       continue
c     # negate the normal velocity:
      do 336 j=1,nyb
         do 336 i=1,nrow
            val(3,i,j) = -val(3,i,j)
  336    continue
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      if (yhi_patch .le. yupper+hymarg) then
c        # not a physical boundary --  no cells at this edge lies
c        # outside the physical bndry.
c        # values are set elsewhere in amr code.
         go to 499
         endif
c
c     # number of grid cells lying outside physical domain:
      nyt = (yhi_patch - yupper + hymarg)/hy
      jbeg = max0(ncol-nyt+1, 1)
c
      go to (400,410,420,430) mthbc(4)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) 
     &   '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2amr'
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do 415 j=jbeg,ncol
         do 415 i=1,nrow
            do 415 m=1,meqn
               val(m,i,j) =  val(m,i,jbeg-1)
  415       continue
      go to 499

  420 continue
c     # periodic:   handled elsewhere in amr
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 j=jbeg,ncol
         do 435 i=1,nrow
            do 435 m=1,meqn
               val(m,i,j) =  val(m,i,2*jbeg-1-j)
  435       continue
c     # negate the normal velocity:
      do 436 j=jbeg,ncol
         do 436 i=1,nrow
            val(3,i,j) = -val(3,i,j)
  436    continue
      go to 499

  499 continue

      return
      end

