c
!> For a rectangle defined on level **level** and bound by 
!! **ilo**, **ihi**, **jlo**, **jhi**, find intersecting grids at the same level
!! and copy data from intersecting grids to solution and auxiliary variables 
!! on the rectangle, stored in **val** and **aux** arrays.
c --------------------------------------------------------------------
c
       subroutine icall(val,aux,nrow,ncol,nvar,naux,
     .                  ilo,ihi,jlo,jhi,level,iputst,jputst)

       use amr_module
       implicit double precision (a-h, o-z)

       dimension val(nvar,nrow,ncol)
       dimension aux(naux,nrow,ncol)

       logical sticksout


c NEW INDEX ORDERING
       iadd   (ivar,i,j) = loc    + ivar-1 + nvar*((j-1)*mitot+i-1)
       iaddaux(ivar,i,j) = locaux + ivar-1 + naux*((j-1)*mitot+i-1)

c ::::::::::::::::::::::::::: ICALL :::::::::::::::::::::::::::::::
c
c    find intersecting grids at the same level. copy data from
c    intersecting grids to both val and aux arrays.
c
c    use larger definition of grids here - boundary data already in.
c    aux arrays also enlarged size.
c
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
c$$$          ixlo = max(iglo-nghost,ilo)
c$$$          ixhi = min(ighi+nghost,ihi)
c$$$          jxlo = max(jglo-nghost,jlo)
c$$$          jxhi = min(jghi+nghost,jhi)
c  how did ghost cells get in the allowable region? They are not filled
c  (since we may be interpolating from newly filled grids, not just grids
c  that have been primed with bcs to be advanced.
          ixlo = max(iglo,ilo)
          ixhi = min(ighi,ihi)
          jxlo = max(jglo,jlo)
          jxhi = min(jghi,jhi)


          if (ixlo .le. ixhi .and. jxlo .le. jxhi) then
              loc  = node(store1,mptr)
              locaux = node(storeaux,mptr)
              nx   = ighi - iglo + 1
              ny   = jghi - jglo + 1
              mitot = nx + 2*nghost
              mjtot = ny + 2*nghost
              do 30 j    = jxlo, jxhi
              do 30 i    = ixlo, ixhi
              do 20 ivar = 1, nvar
                  ialloc  =  iadd(ivar,i-iglo+nghost+1,j-jglo+nghost+1)
                  val(ivar,i-ilo+iputst,j-jlo+jputst)  =  alloc(ialloc)
 20           continue
              do 25 iaux = 1, naux
                  ialloc = iaddaux(iaux,i-iglo+nghost+1,j-jglo+nghost+1)
                  aux(iaux,i-ilo+iputst,j-jlo+jputst)  =  alloc(ialloc)
 25           continue
 30           continue
          endif
          mptr = node(levelptr, mptr)
          go to 10

 99   continue

c   if cells stick out of domain but not periodic then set elsewhere 
c   either setaux  and bc2amr. (called from routine that called this, e.g.
c   saveqc or filval)

        return
        end
