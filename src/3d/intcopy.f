c
c --------------------------------------------------------------------
c
       subroutine intcopy(val,mitot,mjtot,mktot,nvar,
     &                    ilo,ihi,jlo,jhi,klo,khi,level,
     &                    iputst,jputst,kputst)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension val(nvar,mitot,mjtot,mktot)


      iadd(ivar,i,j,k)   = loc    +    (ivar-1)
     &                            +    (i-1)*nvar
     &                            +    (j-1)*nvar*mi
     &                            +    (k-1)*nvar*mi*mj

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
          kglo = node(ndklo,mptr)
          kghi = node(ndkhi,mptr)

c         # does it intersect?
          ixlo = max(iglo,ilo)
          ixhi = min(ighi,ihi)
          jxlo = max(jglo,jlo)
          jxhi = min(jghi,jhi)
          kxlo = max(kglo,klo)
          kxhi = min(kghi,khi)

          if ((ixlo .le. ixhi .and. jxlo .le. jxhi) .and.
     c                             (kxlo .le. kxhi)) then
              loc  = node(store1,mptr)
              nx   = ighi - iglo + 1
              ny   = jghi - jglo + 1
              nz   = kghi - kglo + 1
              mi   = nx + 2*nghost
              mj   = ny + 2*nghost
              mk   = nz + 2*nghost
              do 20 k    = kxlo, kxhi
              do 20 j    = jxlo, jxhi
              do 20 i    = ixlo, ixhi
              do 40 ivar = 1, nvar
                  val(ivar,iputst+i-ilo,jputst+j-jlo,kputst+k-klo) =
     1                alloc(iadd(ivar,i-iglo+nghost+1,j-jglo+nghost+1,
     2                                           k-kglo+nghost+1))
 40           continue
 20           continue
          endif
          mptr = node(levelptr, mptr)
          go to 10

 99   return
      end
