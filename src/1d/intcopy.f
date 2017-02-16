c
c --------------------------------------------------------------------
c
       subroutine intcopy(val,mitot,nvar,ilo,ihi,level,iputst)

      use amr_module
       implicit double precision (a-h, o-z)

       dimension val(nvar,mitot)


c   OLD INDEXING
c      iadd(i,ivar) = loc + i - 1
c   NEW INDEXING ORDER SWITCHED
       iadd(ivar,i) = loc + ivar-1 + nvar*(i-1)

c ::::::::::::::::::::::::::: INTCOPY :::::::::::::::::::::::::::::::
c
c    find intersecting grids at the same level. copy data from
c    old grid to val. 
c    old grid has "nghost" ghost cells - passed in nodal common block.
c    new grid has no ghost cells - indices describe entire patch.
c    iputst: where to copy values into. may not be in
c                    location corresponding to ilo,ihi if
c                    the patch has been periodically wrapped.
c    
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


       mptr = lstart(level)

 10    if (mptr .eq. 0) go to 99
          iglo = node(ndilo,mptr)
          ighi = node(ndihi,mptr)

c         # does it intersect?
          ixlo = max(iglo,ilo)
          ixhi = min(ighi,ihi)

          if (ixlo .le. ixhi) then
              loc  = node(store1,mptr)
              nx   = ighi - iglo + 1
              mi   = nx + 2*nghost
              do 20 ivar = 1, nvar
              do 30 i    = ixlo, ixhi
                 val(ivar,iputst+i-ilo) =
     1               alloc(iadd(ivar,i-iglo+nghost+1))
 30           continue
 20           continue
          endif
          mptr = node(levelptr, mptr)
          go to 10

 99   return
      end
