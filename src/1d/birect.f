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
c
10    continue
      cxlo    = rnode(cornxlo,mptr)
      cxhi    = rnode(cornxhi,mptr)
      nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
      minsize = 2*nghost
c
c check number of rows - if too many, bisect grid with vertical
c line down the middle. make sure new grid corners are anchored
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

        rnode(cornxlo,mptrnx)    = cxmid
        rnode(cornxhi,mptrnx)    = cxhi
        rnode(timemult,mptrnx)   = rnode(timemult,mptr)
        node(nestlevel,mptrnx)   = node(nestlevel,mptr)

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
