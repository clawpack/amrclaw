c
c --------------------------------------------------------------------
c
       subroutine icall(val,aux,nrow,ncol,nfil,nvar,naux,
     .                  ilo,ihi,jlo,jhi,klo,khi,level,
     .                  iputst,jputst,kputst)

      use amr_module
      implicit double precision (a-h, o-z)

      dimension val(nvar,nrow,ncol,nfil)
      dimension aux(naux,nrow,ncol,nfil)


       iadd   (ivar,i,j,k) = loc    +    (ivar-1)
     &                              +    (i-1)*nvar
     &                              +    (j-1)*nvar*mitot
     &                              +    (k-1)*nvar*mitot*mjtot
       iaddaux(ivar,i,j,k) = locaux +    (ivar-1)
     &                              +    (i-1)*naux
     &                              +    (j-1)*naux*mitot
     &                              +    (k-1)*naux*mitot*mjtot

c ::::::::::::::::::::::::::: ICALL :::::::::::::::::::::::::::::::
c
c    find intersecting grids at the same level. copy data from
c    intersecting grids to both val and aux arrays.
c
c    use larger definition of grids here - boundary data already in.
c    aux arrays also enlarged size.
c
c    iputst, jputst, kputst: where to copy values into. may not be in
c                            location corresponding to ilo,ihi,etc. if
c                            the patch has been periodically wrapped.

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
!--         ixlo = max(iglo-nghost,ilo)
!--         ixhi = min(ighi+nghost,ihi)
!--         jxlo = max(jglo-nghost,jlo)
!--         jxhi = min(jghi+nghost,jhi)
!--         kxlo = max(kglo-nghost,klo)
!--         kxhi = min(kghi+nghost,khi)
c  how did ghost cells get in the allowable region? They are not filled
c  (since we may be interpolating from newly filled grids, not just grids
c  that have been primed with bcs to be advanced.
         ixlo = max(iglo,ilo)
         ixhi = min(ighi,ihi)
         jxlo = max(jglo,jlo)
         jxhi = min(jghi,jhi)
         kxlo = max(kglo,klo)
         kxhi = min(kghi,khi)

         if ((ixlo .le. ixhi .and. jxlo .le. jxhi) .and.
     &       (                     kxlo .le. kxhi)) then
           loc  = node(store1,mptr)
           locaux = node(storeaux,mptr)
           nx   = ighi - iglo + 1
           ny   = jghi - jglo + 1
           nz   = kghi - kglo + 1
           mitot = nx + 2*nghost
           mjtot = ny + 2*nghost
           mktot = nz + 2*nghost
           do 30 k    = kxlo, kxhi
           do 30 j    = jxlo, jxhi
           do 30 i    = ixlo, ixhi
           do 20 ivar = 1, nvar
              ialloc  =  iadd(ivar,i-iglo+nghost+1,j-jglo+nghost+1,
     &                                             k-kglo+nghost+1)
              val(ivar,i-ilo+iputst,j-jlo+jputst,
     &                              k-klo+kputst) = alloc(ialloc)
 20     continue
           do 25 iaux = 1, naux
             ialloc = iaddaux(iaux,i-iglo+nghost+1,j-jglo+nghost+1,
     &                                             k-kglo+nghost+1)
             aux(iaux,i-ilo+iputst,j-jlo+jputst,
     &                             k-klo+kputst) = alloc(ialloc)
 25        continue
 30       continue
       endif
       mptr = node(levelptr, mptr)
       go to 10

 99   return
      end
