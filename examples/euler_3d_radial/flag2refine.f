c
c -------------------------------------------------------------------
      subroutine flag2refine(mx,my,mz,mbc,meqn,maux,xlower,ylower,
     &      zlower,dx,dy,dz,t,
     &      level,tolsp,q,aux,amrflags,DONTFLAG,DOFLAG)
c     ------------------------------------------------------------------
c     -

c
c     ::::::::::::::::::::: flag2refine
c     ::::::::::::::::::::::::::::::::::
c
c     User routine to control flagging of points for refinement.
c
c     Default version computes spatial difference dq in each direction
c     and
c     for each component of q and flags any point where this is greater
c     than
c     the tolerance tolsp.  This is consistent with what the routine
c     errsp did in
c     earlier versions of amrclaw (4.2 and before).
c
c     This routine can be copied to an application directory and
c     modified to
c     implement some other desired refinement criterion.
c
c     The logical function allowflag(x,y,t,level) is called to check
c     whether
c     further refinement at this level is allowed at this particular
c     location
c     and time.  The default library version of this routine returns
c     .true.
c     for all arguments.  Copy that routine to the application directory
c     and
c     modify it if needed to restrict the region where refinement is
c     allowed.
c
c     Points may also be flagged for refining based on a Richardson
c     estimate
c     of the error, obtained by comparing solutions on the current grid
c     and a
c     coarsened grid.  Points are flagged if the estimated error is
c     larger than
c     the parameter tol in amr2ez.data, provided tol>0.  If tol<=0 then
c     the coarsening and Richardson estimation is not performed!
c     This is a change from previous versions (4.2 and before) of
c     amrclaw.
c     Note: in previous versions, the routine errf1 used a function
c     allowed(x,y,level) that has been replaced by the allowflag.  This
c     new
c     function is also used in Richardson estimation if that is invoked.
c
c
c     q   = grid values including ghost cells (bndry vals at specified
c     time have already been set, so can use ghost cell values too)
c
c     aux   = aux array on this grid patch
c
c     amrflags  = array to be flagged with either the value
c     DONTFLAG (no refinement needed)  or
c     DOFLAG   (refinement desired)
c
c     tolsp = tolerance specified by user in input file amr2ez.data,
c     used in default
c     version of this routine as a tolerance for spatial differences.

c     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit none

      integer mx,my,mz,mbc,meqn,maux
      double precision xlower, ylower, zlower, dx,dy,dz,t, tolsp
      double precision DONTFLAG, DOFLAG
      integer level
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      double precision amrflags(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
      logical  allowflag
      external allowflag

      integer m,i,j,k
      double precision x,y,z
      double precision dp(3), dpi, dpj, dpk, pf,dpmax


c     # loop over interior points on this grid:
      do 30 k = 1,mz
         z = zlower + (k-0.5d0)*dz

         do 20 j = 1,my
            y = ylower + (j-0.5d0)*dy
            do 10 i = 1,mx
               x = xlower + (i-0.5d0)*dx

               amrflags(i,j,k) = DONTFLAG

c               if (allowflag(x,y,z,t,level)) then

c                 # check to see if we should flag this point for
c                 refinement.
c                 # Here the default test is taken from errsp.f in
c                 previous
c                 # versions of amrclaw -- flag this point if dq >
c                 tolsp:

                 dp(1) = pf(q(1,i+1,j,k)) - pf(q(1,i-1,j,k))
                 dp(2) = pf(q(1,i,j+1,k)) - pf(q(1,i,j-1,k))
                 dp(3) = pf(q(1,i,j,k+1)) - pf(q(1,i,j,k-1))

                 dpmax = 0.0
                 do m = 1,3
                    dpmax  = max(abs(dp(m)), dpmax)
                 enddo

                 if (abs(dpmax) .gt. tolsp) then
                    amrflags(i,j,k) = DOFLAG
                 endif

c               endif

   10          continue
   20       continue
   30    continue

         return
         end


      double precision function pf(q)
      implicit none

      double precision q(5)
      double precision press, mom2, e, rho

      double precision gamma, gamma1
      common /cparam/  gamma, gamma1

      rho = q(1)
      mom2 = q(2)**2  + q(3)**2 + q(4)**2
      pf = gamma1*(q(5) - 0.5d0*mom2/rho)


      end
