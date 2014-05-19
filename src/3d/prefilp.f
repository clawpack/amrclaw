c
c --------------------------------------------------------------
c
      recursive subroutine prefilrecur(level,nvar,
     1                                 valbig,aux,naux,time,
     2                                 mitot,mjtot,mktot,
     3                                 nrowst,ncolst,nfilst,
     4                                 ilo,ihi,jlo,jhi,klo,khi)

      use amr_module
      implicit double precision (a-h,o-z)


      dimension valbig(nvar,mitot,mjtot,mktot)
      dimension    aux(naux,mitot,mjtot,mktot)
      dimension ist(3), iend(3), ishift(3)
      dimension jst(3), jend(3), jshift(3)
      dimension kst(3), kend(3), kshift(3)
  

c     dimension scratch(max(mitot,mjtot,mktot)*nghost*nvar)
c     dimension scratchaux(max(mitot,mjtot,mktot)*nghost*naux)

c      iadd(ivar,i,j,k) = locflip + ivar - 1 + nvar*(i-1) + nvar*nr*(j-1)
c     &                           + nvar*nr*nc*(k-1)
c      iaddscratch(ivar,i,j,k) = ivar + nvar*(i-1) + nvar*nr*(j-1)
c     &                               + nvar*nr*nc*(k-1)

c
c  :::::::::::::: PREFILRECUR :::::::::::::::::::::::::::::::::::::::::::
c     For periodic boundary conditions more work needed to fill the
c     piece of the boundary. This routine was
c     called because the patch sticks out of the domain,
c     and has periodic bc.s preprocess the patch before calling
c     filrecur to shift the patch periodically back into the domain.
c
c     Inputs to this routine:
c     xl, xr, yb, yt, zf, zr = the location in physical space of
c     corners of a patch.
c     ilo,ihi,jlo,jhi,klo,khi = the location in index space of this patch.
c
c     Outputs from this routine:
c     The values around the border of the grid are inserted
c     directly into the enlarged valbig array for this piece.
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     # will divide patch into 27 possibilities (some empty): 
c       x sticks out left, x interior, x sticks out right
c       same for y and z. for example, the max. would be
c       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)

      if (xperdom) then
         ist(1)  = ilo
         ist(2)  = 0
         ist(3)  = iregsz(level)
         iend(1) = -1
         iend(2) = iregsz(level)-1
         iend(3) = ihi
         ishift(1) = iregsz(level)
         ishift(2) = 0
         ishift(3) = -iregsz(level)
      else
         ist(1)    = iregsz(level)
         ist(2)    = ilo
         ist(3)    = iregsz(level)
         iend(1)   = -1
         iend(2)   = ihi
         iend(3)   = -1
         ishift(1) = 0
         ishift(2) = 0
         ishift(3) = 0
      endif

      if (yperdom) then
         jst(1)  = jlo
         jst(2)  = 0
         jst(3)  = jregsz(level)
         jend(1) = -1
         jend(2) = jregsz(level)-1
         jend(3) = jhi
         jshift(1) = jregsz(level)
         jshift(2) = 0
         jshift(3) = -jregsz(level)
      else
         jst(1)    = jregsz(level)
         jst(2)    = jlo
         jst(3)    = jregsz(level)
         jend(1)   = -1
         jend(2)   = jhi
         jend(3)   = -1
         jshift(1) = 0
         jshift(2) = 0
         jshift(3) = 0
      endif

      if (zperdom) then
         kst(1)  = klo
         kst(2)  = 0
         kst(3)  = kregsz(level)
         kend(1) = -1
         kend(2) = kregsz(level)-1
         kend(3) = khi
         kshift(1) = kregsz(level)
         kshift(2) = 0
         kshift(3) = -kregsz(level)
      else
         kst(1)    = kregsz(level)
         kst(2)    = klo
         kst(3)    = kregsz(level)
         kend(1)   = -1
         kend(2)   = khi
         kend(3)   = -1
         kshift(1) = 0
         kshift(2) = 0
         kshift(3) = 0
      endif


        do 30 i = 1, 3
         i1 = max(ilo,  ist(i))
         i2 = min(ihi, iend(i))
         if (i1 .gt. i2) go to 30

         do 20 j = 1, 3
           j1 = max(jlo,  jst(j))
           j2 = min(jhi, jend(j))
           if (j1 .gt. j2) go to 20

           do 10 k = 1, 3
             k1 = max(klo,  kst(k))
             k2 = min(khi, kend(k))

           if (k1 <= k2) then ! part of patch in this region
              iputst = (i1 - ilo) + nrowst
              jputst = (j1 - jlo) + ncolst
              kputst = (k1 - klo) + nfilst

              kuse1 = k1+kshift(k)
              kuse2 = k2+kshift(k)

              call filrecur(level,nvar,valbig,aux,naux,time,
     1                      mitot,mjtot,mktot,
     2                      iputst,jputst,kputst,
     3                      i1+ishift(i),i2+ishift(i),
     4                      j1+jshift(j),j2+jshift(j),
     5                      kuse1,kuse2)               

           end if

 10      end do
 20     end do
 30    end do

       return
       end
