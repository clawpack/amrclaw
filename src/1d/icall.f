c
c --------------------------------------------------------------------
c
       subroutine icall(val,aux,nrow,nvar,naux,
     .                  ilo,ihi,level,iputst)

       use amr_module
       implicit double precision (a-h, o-z)

       dimension val(nvar,nrow)
       dimension aux(naux,nrow)


c NEW INDEX ORDERING
       iadd   (ivar,i) = loc    + ivar-1 + nvar*(i-1)
       iaddaux(ivar,i) = locaux + ivar-1 + naux*(i-1)

c ::::::::::::::::::::::::::: ICALL :::::::::::::::::::::::::::::::
c
c    find intersecting grids at the same level. copy data from
c    intersecting grids to both val and aux arrays.
c
c    use larger definition of grids here - boundary data already in.
c    aux arrays also enlarged size.
c
c    iputst: where to copy values into. may not be in
c            location corresponding to ilo,ihi if
c            the patch has been periodically wrapped.

c    
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


       mptr = lstart(level)

 10    if (mptr .eq. 0) go to 99
          iglo = node(ndilo,mptr) 
          ighi = node(ndihi,mptr)

c         # does it intersect?
c$$$          ixlo = max(iglo-nghost,ilo)
c$$$          ixhi = min(ighi+nghost,ihi)
c  how did ghost cells get in the allowable region? They are not filled
c  (since we may be interpolating from newly filled grids, not just grids
c  that have been primed with bcs to be advanced.
          ixlo = max(iglo,ilo)
          ixhi = min(ighi,ihi)


          if (ixlo .le. ixhi) then
              loc  = node(store1,mptr)
              locaux = node(storeaux,mptr)
              nx   = ighi - iglo + 1
              mitot = nx + 2*nghost
              do 30 i    = ixlo, ixhi
              do 20 ivar = 1, nvar
                  ialloc  =  iadd(ivar,i-iglo+nghost+1)
                  val(ivar,i-ilo+iputst)  =  alloc(ialloc)
 20           continue
              do 25 iaux = 1, naux
                  ialloc = iaddaux(iaux,i-iglo+nghost+1)
                  aux(iaux,i-ilo+iputst)  =  alloc(ialloc)
 25           continue
 30           continue
          endif
          mptr = node(levelptr, mptr)
          go to 10

 99   continue

c   if cells stick out of domain but not periodic then set elsewhere 
c   either setaux  and bc1amr. (called from routine that called this, e.g.
c   saveqc or filval)

        return
        end
