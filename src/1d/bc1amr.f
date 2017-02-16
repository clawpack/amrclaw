
c
c ------------------------------------------------------------------
c
      subroutine bc1amr(val,aux,nrow,meqn,naux,
     1                  hx, level, time,
     2                  xlo_patch, xhi_patch)
 
c
c
c :::::::::: bc1amr ::::::::::::::::::::::::::::::::::::::::::::::;
c
c     Take a grid patch with mesh width hx of dimensions nrow
c     and set the values of any piece of
c     of the patch which extends outside the physical domain 
c     using the boundary conditions. 
c
c     ------------------------------------------------
c     # Standard boundary condition choices for amr in clawpack
c
c     # At each boundary  k = 1 (left),  2 (right):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary conditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) component of q.
c     #            =  5  sphere bcs (left half maps to right half of same 
c     #                  side, and vice versa), as if domain folded in half
c     ------------------------------------------------
c
c     The edges of the grid patch are at
c        xlo_patch  -- left edge
c        xhi_patch --  right edge
c
c     The physical domain itself is a rectangle bounded by
c        xlower  -- left edge
c        xupper  -- right edge
c     
c     the picture is the following: 
c            (xlower)                    (xupper)
c        |______|_________|_________________|
c        |      |         |                 |
c   (xlo_patch)       (xhi_patch)
c        
c
c     Any cells that lie outside the physical domain are ghost cells whose
c     values should be set in this routine.  This is tested for by comparing
c     xlo_patch with xlower to see if values need to be set at the left, as in
c     the figure above, and similarly at the other boundary.
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

      use amr_module, only: mthbc,xlower,xupper

      implicit none

      integer nrow,meqn,naux,level
      real(kind=8) :: hx,time, hxmarg
      real(kind=8) :: xlo_patch,xhi_patch
      integer nxl,nxr,ibeg,i,m
      real(kind=8) :: val(meqn,nrow), aux(naux,nrow)
      

      hxmarg = hx*.01

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
         do 115 i=1,nxl
            do 115 m=1,meqn
               val(m,i) = val(m,nxl+1)
  115       continue
      go to 199

  120 continue
c     # periodic:   handled elsewhere in amr
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
         do 135 i=1,nxl
            do 135 m=1,meqn
               val(m,i) = val(m,2*nxl+1-i)
  135       continue
c     # negate the normal velocity:
         do 136 i=1,nxl
            val(2,i) = -val(2,i)
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
         do 215 i=ibeg,nrow
            do 215 m=1,meqn
               val(m,i) = val(m,ibeg-1)
  215       continue
      go to 299
  
  220 continue
c     # periodic:   handled elsewhere in amr
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
         do 235 i=ibeg,nrow
            do 235 m=1,meqn
               val(m,i) = val(m,2*ibeg-1-i)
  235       continue
c     # negate the normal velocity:
         do 236 i=ibeg,nrow
            val(2,i) = -val(2,i)
  236    continue
      go to 299
   
  299 continue

      return
      end

