c
!> For a rectangle that is on level **level**, described by
!! **ilo**, **ihi**, **jlo**, **jhi** and made up by 
!! **mitot** \f$ \times\f$ **mjtot** cells, copy solution from 
!! all (OLD) grids on level **level** to **val**, which stores
!! solution on this rectangle.
!! Some portions of the rectangle may not be filled since they do not
!! overlap with any level **level** grids.
!!  
c --------------------------------------------------------------------
c
       subroutine intcopy(val,mitot,mjtot,nvar,ilo,ihi,jlo,jhi,level,
     &                    iputst,jputst)

      use amr_module
       implicit double precision (a-h, o-z)

       dimension val(nvar,mitot,mjtot)


c   OLD INDEXING
c      iadd(i,j,ivar) = loc + i - 1 + mi*((ivar-1)*mj+j-1)
c   NEW INDEXING ORDER SWITCHED
       iadd(ivar,i,j) = loc + ivar-1 + nvar*((j-1)*mi+i-1)

c ::::::::::::::::::::::::::: INTCOPY :::::::::::::::::::::::::::::::
c
c    find intersecting grids at the same level. copy data from
c    old grid to val. 
c    old grid has "nghost" ghost cells - passed in nodal common block.
c    new grid has no ghost cells - indices describe entire patch.
c    iputst, jputst: where to copy values into. may not be in
c                    location corresponding to ilo,ihi,etc. if
c                    the patch has been periodically wrapped.
c    
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


       mptr = lstart(level)

 10    if (mptr .eq. 0) go to 99
          iglo = node(ndilo,mptr)
          ighi = node(ndihi,mptr)
          jglo = node(ndjlo,mptr)
          jghi = node(ndjhi,mptr)

c         # does it intersect?
          ixlo = max(iglo,ilo)
          ixhi = min(ighi,ihi)
          jxlo = max(jglo,jlo)
          jxhi = min(jghi,jhi)

          if (ixlo .le. ixhi .and. jxlo .le. jxhi) then
              loc  = node(store1,mptr)
              nx   = ighi - iglo + 1
              ny   = jghi - jglo + 1
              mi   = nx + 2*nghost
              mj   = ny + 2*nghost
              do 20 j    = jxlo, jxhi
              do 21 ivar = 1, nvar
              do 30 i    = ixlo, ixhi
                 val(ivar,iputst+i-ilo,jputst+j-jlo) =
     1               alloc(iadd(ivar,i-iglo+nghost+1,j-jglo+nghost+1))
 30           continue
 21           continue
 20           continue
          endif
          mptr = node(levelptr, mptr)
          go to 10

 99   return
      end
