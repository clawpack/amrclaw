c
c ---------------------------------------------------------
c
      subroutine outmsh(mptr,outgrd,nvar,naux)

      use amr_module
      implicit double precision (a-h,o-z)

      logical  outgrd
c
c ::::::::::::::::::::: OUTMSH :::::::::::::::::::::::::::::::::::::
c
c outmsh - output the grid descriptor and optionally the values on
c          the grid (for a single grid - see "outtre" for outputing
c          a subtree)
c input parameters:
c    mptr   - ptr to grid descriptor
c    outgrd - if true, output value on grid
c special case
c     if grid has level < 1, nothing is printed. (forest pointer
c has level = 0).  this simplifies logic of outtre; also, any grid
c with level <= 0 is not properly set (the forest pointer
c is used only to provide pointers into the tree of coarsest grids)
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
100   format(1x,47h+----------------------------------------------,
     *30h-----------------------------+)
c
      lev = node(nestlevel,mptr)
c
      write(outunit,100)
      write(outunit,101) mptr
101   format(1x,10h! grid no:,i6,60x,1h!)
      write(outunit,102) node(nestlevel,mptr),rnode(timemult,mptr),
     .             node(levelptr,mptr)
102   format(1x,1h!,11h nestlevel=,i3,12h, time mult=,f8.5,
     1       13h, level ptr =,i6,22x,1h!)
      write(outunit,103) node(store1,mptr),node(store2,mptr),
     1                   node(cfluxptr,mptr),node(ffluxptr,mptr)
 103  format(1x,'! storage locs =',2i12,'  bndry locs =',2i11,1h!)
      write(outunit,104)
      write(outunit,111) rnode(cornxhi,mptr),
     2                   rnode(cornyhi,mptr),
     3                   rnode(cornzhi,mptr)
      write(outunit,111) rnode(cornxlo,mptr),
     2                   rnode(cornylo,mptr),
     3                   rnode(cornzlo,mptr)
      write(outunit,112)
      write(outunit,113) node(ndihi,mptr),
     2                   node(ndjhi,mptr),
     3                   node(ndkhi,mptr)
      write(outunit,113) node(ndilo,mptr),
     2                   node(ndjlo,mptr),
     3                   node(ndklo,mptr)
112   format(1x,23h! integer index space :,53x,1h!)
113   format(1x,2h! ,18x,1h(,i8,2h, ,i8,2h, ,i8,1h),26x,1h!)
104   format(1x,40h! corners of rectangular parallelopiped:,36x,1h!)
111   format(1x,2h! ,18x,1h(,f8.5,2h, ,f8.5,2h, ,f8.5,1h),26x,1h!)
      write(outunit,105) hxposs(lev),hyposs(lev),hzposs(lev),possk(lev)
105   format(1x,7h! hrow=,f9.6,7h, hcol=,f9.6,7h, hfil=,f9.6,
     .          8h, ktime=,f9.6,11x,1h!)
      write(outunit,100)
c
      if (.not. outgrd) go to 99
c output the grid
      mitot   = node(ndihi,mptr) - node(ndilo,mptr) + 1 + 2*nghost
      mjtot   = node(ndjhi,mptr) - node(ndjlo,mptr) + 1 + 2*nghost
      mktot   = node(ndkhi,mptr) - node(ndklo,mptr) + 1 + 2*nghost
      loc     = node(store1,mptr)
      locaux  = node(storeaux,mptr)
      call outval(alloc(loc),nvar,mitot,mjtot,mktot,mptr,outgrd,
     1            naux,alloc(locaux))
 99   return
      end
