c
c --------------------------------------------------------------
c
      recursive subroutine prefilrecur(level,nvar,
     1                                 valbig,aux,naux,time,
     2                                 mitot,mjtot,mktot,
     3                                 nrowst,ncolst,nfilst
     4                                 ilo,ihi,jlo,jhi,klo,khi)

      implicit double precision (a-h,o-z)

      include  "call.i"

      dimension valbig(mitot,mjtot,mktot,nvar)
      dimension    aux(mitot,mjtot,mktot,naux)
      dimension ist(3), iend(3), jst(3), jend(3), kst(3), kend(3)
      dimension ishift(3), jshift(3), kshift(3)

      dimension scratch(max(mitot,mjtot,mktot)*nghost*nvar)
      dimension scratchaux(max(mitot,mjtot,mktot)*nghost*naux)

      iadd(i,j,k,ivar)  = locflip + i - 1 + nr*(j-1) + nc*nr*(k-1)
     &                            + nf*nc*nr*(ivar-1)
      iaddscratch(i,j,k,ivar)  = i + nr*(j-1) + nc*nr*(k-1)
     &                             + nf*nc*nr*(ivar-1)

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

        ist(1)  = ilo
        ist(2)  = 0
        ist(3)  = iregsz(level)
        iend(1) = -1
        iend(2) = iregsz(level)-1
        iend(3) = ihi

        jst(1)  = jlo
        jst(2)  = 0
        jst(3)  = jregsz(level)
        jend(1) = -1
        jend(2) = jregsz(level)-1
        jend(3) = jhi

        kst(1)  = klo
        kst(2)  = 0
        kst(3)  = kregsz(level)
        kend(1) = -1
        kend(2) = kregsz(level)-1
        kend(3) = khi

        ishift(1) = iregsz(level)
        ishift(2) = 0
        ishift(3) = -iregsz(level)
        jshift(1) = jregsz(level)
        jshift(2) = 0
        jshift(3) = -jregsz(level)
        kshift(1) = kregsz(level)
        kshift(2) = 0
        kshift(3) = -kregsz(level)


        do i = 1, 3
           i1 = max(ilo,  ist(i))
           i2 = min(ihi, iend(i))
           do j = 1, 3
           j1 = max(jlo,  jst(j))
           j2 = min(jhi, jend(j))
           do k = 1, 3
           k1 = max(klo,  kst(j))
           k2 = min(khi, kend(j))

           if ((i1 <= i2) .and. (j1 <= j2) .and. (k1 <= k2)) then ! part of patch in this region
              iputst = (i1 - ilo) + nrowst
              jputst = (j1 - jlo) + ncolst
              kputst = (k1 - klo) + nfilst
              call filrecur(level,nvar,valbig,aux,naux,time,
     1                      mitot,mjtot,mktot,
     2                      iputst,jputst,kputst,
     3                      i1+ishift(i),i2+ishift(i),
     4                      j1+jshift(j),j2+jshift(j),
     5                      k1+kshift(k),k2+kshift(k))

           end if

           end do
           end do
        end do

        return
        end
