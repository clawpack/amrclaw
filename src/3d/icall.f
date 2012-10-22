c
c --------------------------------------------------------------------
c
       subroutine icall(val,aux,nrow,ncol,nfil,nvar,naux,
     .                  ilo,ihi,jlo,jhi,klo,khi,level,
     .                  iputst,jputst,kputst)

       implicit double precision (a-h, o-z)

       dimension val(nrow,ncol,nfil,nvar)
       dimension aux(nrow,ncol,nfil,naux)

       include "call.i"

c2D    iadd   (i,j,ivar) = loc    + i - 1 + mitot*((ivar-1)*mjtot+j-1)
c2D    iaddaux(i,j,ivar) = locaux + i - 1 + mitot*((ivar-1)*mjtot+j-1)
       iadd   (i,j,k,ivar) = loc    +    (i-1)
     &                              +    (j-1)*mitot
     &                              +    (k-1)*mitot*mjtot
     &                              + (ivar-1)*mitot*mjtot*mktot
       iaddaux(i,j,k,ivar) = locaux +    (i-1)
     &                              +    (j-1)*mitot
     &                              +    (k-1)*mitot*mjtot
     &                              + (ivar-1)*mitot*mjtot*mktot

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
	  ixlo = max(iglo-nghost,ilo)
	  ixhi = min(ighi+nghost,ihi)
	  jxlo = max(jglo-nghost,jlo)
	  jxhi = min(jghi+nghost,jhi)
	  kxlo = max(kglo-nghost,klo)
	  kxhi = min(kghi+nghost,khi)

	  if ((ixlo .le. ixhi .and. jxlo .le. jxhi) .and.
     &        (                     kxlo .le. kxhi)) then
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
		  ialloc  =  iadd(i-iglo+nghost+1,j-jglo+nghost+1,
     &                                            k-kglo+nghost+1,ivar)
		  val(i-ilo+iputst,j-jlo+jputst,
     &                             k-klo+kputst,ivar)  =  alloc(ialloc)
 20           continue
              do 25 iaux = 1, naux
                  ialloc = iaddaux(i-iglo+nghost+1,j-jglo+nghost+1,
     &                                             k-kglo+nghost+1,iaux)
		  aux(i-ilo+iputst,j-jlo+jputst,
     &                             k-klo+kputst,iaux)  =  alloc(ialloc)
 25           continue
 30           continue
	  endif
	  mptr = node(levelptr, mptr)
	  go to 10

 99   return
      end
