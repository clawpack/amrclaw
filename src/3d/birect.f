c
c --------------------------------------------------
c
      subroutine birect(mptr1)
c
      use amr_module
      implicit double precision (a-h,o-z)

c
c :::::::::::::  BIRECT :::::::::::::::::::::::::::::::::::::::
c check each grid, starting with mptr1 (either newstl or lstart)
c to see that it has no more than max1d points in either dimensions.
c needed so that scratch array space in stepgrid not exceeded.
c
c also check for too small grids - but has never happened.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      mptr  = mptr1
      level = node(nestlevel,mptr)
      hx    = hxposs(level)
      hy    = hyposs(level)
      hz    = hzposs(level)
c
10    continue
      cxlo    = rnode(cornxlo,mptr)
      cxhi    = rnode(cornxhi,mptr)
      cylo    = rnode(cornylo,mptr)
      cyhi    = rnode(cornyhi,mptr)
      czlo    = rnode(cornzlo,mptr)
      czhi    = rnode(cornzhi,mptr)

      nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
      ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
      nz      = node(ndkhi,mptr) - node(ndklo,mptr) + 1
      minsize = 2*nghost
C  #      if ((nx .lt. minsize .or. ny .lt. minsize .or. nz .lt. minsize)
C  #     .                   .and. (level .lt. lfine)) then
C  #         write(6,*) " made too small a grid - found in birect"
C  #         write(6,*) " fix this problem after all"
C  #         stop
C  #      endif
c
c check number of rows first - if too many, bisect grid with constant x
c plane down the middle. make sure new grid corners are anchored
c on coarse grid point. make sure if bisecting coarse grid that
c new grids have even number of points
c
      if (nx + 2*nghost .gt. max1d) then

         nxl    = nx/2
         if (level .gt. 1) then
           lratio = intratx(level-1)
        else
           lratio = 2
        endif
        nxl = (nxl/lratio)*lratio
        nxr    = nx - nxl
        cxmid  = cxlo + nxl*hx

        mptrnx = nodget()
        node(levelptr,mptrnx) = node(levelptr,mptr)
        node(levelptr,mptr)   = mptrnx

        rnode(cornxhi,mptr) = cxmid

        node(ndihi,mptrnx)  = node(ndihi,mptr)
        node(ndihi,mptr)    = node(ndilo,mptr) + nxl - 1
        node(ndilo,mptrnx)  = node(ndihi,mptr) + 1
        node(ndjhi,mptrnx)  = node(ndjhi,mptr)
        node(ndjlo,mptrnx)  = node(ndjlo,mptr)
        node(ndkhi,mptrnx)  = node(ndkhi,mptr)
        node(ndklo,mptrnx)  = node(ndklo,mptr)

        rnode(cornxlo,mptrnx)    = cxmid
        rnode(cornylo,mptrnx)    = cylo
        rnode(cornzlo,mptrnx)    = czlo
        rnode(cornxhi,mptrnx)    = cxhi
        rnode(cornyhi,mptrnx)    = cyhi
        rnode(cornzhi,mptrnx)    = czhi
        rnode(timemult,mptrnx)   = rnode(timemult,mptr)
        node(nestlevel,mptrnx)   = node(nestlevel,mptr)

        go to 10
c
c check number of columns next - if too many, bisect grid with constant y
c plane down the middle
c
      else if (ny + 2*nghost .gt. max1d) then

         nyl    = ny/2
         if (level .gt. 1) then
            lratio = intraty(level-1)
         else
            lratio = 2
         endif
         nyl = (nyl/lratio)*lratio
         nyr    = ny - nyl
         cymid  =  cylo + nyl*hy

         mptrnx = nodget()
         node(levelptr,mptrnx) = node(levelptr,mptr)
         node(levelptr,mptr)   = mptrnx

         rnode(cornyhi,mptr)   = cymid

         node(ndjhi,mptrnx) = node(ndjhi,mptr)
         node(ndjhi,mptr)   = node(ndjlo,mptr) + nyl - 1
         node(ndjlo,mptrnx) = node(ndjhi,mptr) + 1
         node(ndkhi,mptrnx) = node(ndkhi,mptr)
         node(ndklo,mptrnx) = node(ndklo,mptr)
         node(ndihi,mptrnx) = node(ndihi,mptr)
         node(ndilo,mptrnx) = node(ndilo,mptr)

         rnode(cornxlo,mptrnx)   = cxlo
         rnode(cornylo,mptrnx)   = cymid
         rnode(cornzlo,mptrnx)   = czlo
         rnode(cornxhi,mptrnx)   = cxhi
         rnode(cornyhi,mptrnx)   = cyhi
         rnode(cornzhi,mptrnx)   = czhi
         node(nestlevel,mptrnx)  = node(nestlevel,mptr)
         rnode(timemult,mptrnx)  = rnode(timemult,mptr)
         go to 10
c
c check number of files next - if too many, bisect grid with constant z
c plane down the middle
c
      else if (nz + 2*nghost  .gt. max1d) then

        nzl    = nz/2
        if (level .gt. 1) then
           lratio = intratz(level-1)
        else
           lratio = 2
        endif
        nzl = (nzl/lratio)*lratio
        nzr    = nz - nzl
        czmid  =  czlo + nzl*hz

        mptrnx = nodget()
        node(levelptr,mptrnx) = node(levelptr,mptr)
        node(levelptr,mptr)   = mptrnx

        rnode(cornzhi,mptr)   = czmid

        node(ndkhi,mptrnx) = node(ndkhi,mptr)
        node(ndkhi,mptr)   = node(ndklo,mptr) + nzl - 1
        node(ndklo,mptrnx) = node(ndkhi,mptr) + 1
        node(ndihi,mptrnx) = node(ndihi,mptr)
        node(ndilo,mptrnx) = node(ndilo,mptr)
        node(ndjhi,mptrnx) = node(ndjhi,mptr)
        node(ndjlo,mptrnx) = node(ndjlo,mptr)

        rnode(cornxlo,mptrnx)   = cxlo
        rnode(cornylo,mptrnx)   = cylo
        rnode(cornzlo,mptrnx)   = czmid
        rnode(cornxhi,mptrnx)   = cxhi
        rnode(cornyhi,mptrnx)   = cyhi
        rnode(cornzhi,mptrnx)   = czhi
        node(nestlevel,mptrnx)  = node(nestlevel,mptr)
        rnode(timemult,mptrnx)  = rnode(timemult,mptr)
        go to 10
c
c  grid ok - check the next
c
      else
        mptr = node(levelptr,mptr)
        if (mptr.ne.0) go to 10

      endif

      return
      end
