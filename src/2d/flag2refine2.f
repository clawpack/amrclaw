c
c -------------------------------------------------------------------
      subroutine flag2refine2(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,
     &                       dx,dy,t,level,tolsp,q,aux,
     &                       amrflags,DONTFLAG,DOFLAG)
c -------------------------------------------------------------------

c
c ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
c
c User routine to control flagging of points for refinement.
c
c Default version computes spatial difference dq in each direction and
c for each component of q and flags any point where this is greater than
c the tolerance tolsp.  This is consistent with what the routine errsp did in
c earlier versions of amrclaw (4.2 and before).
c
c This routine can be copied to an application directory and modified to
c implement some other desired refinement criterion.
c
c The logical function allowflag(x,y,t,level) is called to check whether 
c further refinement at this level is allowed at this particular location
c and time.  The default library version of this routine returns .true.
c for all arguments.  Copy that routine to the application directory and
c modify it if needed to restrict the region where refinement is allowed.
c
c Points may also be flagged for refining based on a Richardson estimate
c of the error, obtained by comparing solutions on the current grid and a
c coarsened grid.  Points are flagged if the estimated error is larger than
c the parameter tol in amr2ez.data, provided flag_richardson is .true.,
c otherwise the coarsening and Richardson estimation is not performed!  
c
c This is a change from previous versions (4.2 and before) of amrclaw.
c Note: in previous versions, the routine errf1 used a function
c allowed(x,y,level) that has been replaced by the allowflag.  This new
c function is also used in Richardson estimation if that is invoked.
c
c
c    q   = grid values including ghost cells (bndry vals at specified
c          time have already been set, so can use ghost cell values too)
c
c  aux   = aux array on this grid patch
c
c amrflags  = array to be flagged with either the value
c             DONTFLAG (no refinement needed)  or
c             DOFLAG   (refinement desired)    
c
c tolsp = tolerance specified by user in input file amr2ez.data, used in default
c         version of this routine as a tolerance for spatial differences.

c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit double precision (a-h, o-z)

      dimension   q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension   aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension   amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)
      logical     allowflag
      external    allowflag
      logical     db/.true./
   
c     # loop over interior points on this grid:
      do 20 j = 1,my
        y = ylower + (j-0.5d0)*dy
        do 10 i = 1,mx
          x = xlower + (i-0.5d0)*dx

          amrflags(i,j) = DONTFLAG

          if (allowflag(x,y,t,level)) then
c            # check to see if we should flag this point for refinement.
c            # Here the default test is taken from errsp.f in previous
c            # versions of amrclaw -- flag this point if dq > tolsp:

             dq = 0.d0
             do 5 m = 1,meqn
                dqi = dabs(q(m,i+1,j) - q(m,i-1,j))
                dqj = dabs(q(m,i,j+1) - q(m,i,j-1))
                dq  = dmax1(dq,dqi, dqj)
    5           continue

             if (dq .gt. tolsp) then
                 amrflags(i,j) = DOFLAG
             endif
             endif

 10          continue
 20       continue

      if (db) then
         write(66,*) 'AMR flags from flag2refine:'
         do j = my,1,-1
            write(66,100) (amrflags(i,j),i=1,mx)
 100        format(100f2.0)
         end do   
      endif

      return
      end