c
c --------------------------------------------------------------
c
      recursive subroutine prefilrecur(level,nvar,
     1                                 valbig,auxbig,naux,time,
     2                                 mitot,mjtot,mktot,
     3                                 nrowst,ncolst,nfilst,
     4                                 ilo,ihi,jlo,jhi,klo,khi,
     5                                 iglo,ighi,jglo,jghi,kglo,kghi,
     6                                 patchOnly)

      use amr_module
      implicit double precision (a-h,o-z)


      dimension valbig(nvar,mitot,mjtot,mktot)
      dimension auxbig(naux,mitot,mjtot,mktot)
      dimension ist(3), iend(3), ishift(3)
      dimension jst(3), jend(3), jshift(3)
      dimension kst(3), kend(3), kshift(3)
      logical   patchOnly

    ! dimension scratch patches at largest possible
      dimension  valPatch((ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1)*nvar)  
      dimension  auxPatch((ihi-ilo+1)*(jhi-jlo+1)*(khi-klo+1)*naux)    
      integer    msrc/ -1/  ! indicates not a real grid calling filpatch

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
c     iglo,ighi,jglo,jghi,kglo,kghi = the location in index space of the grid
c                          containing the patch, used to compute where to copy and fill
c     NB: the patch may actually be the whole grid (if called from filpatch for coarser lev)
c
c     Outputs from this routine:
c     The values around the border of the grid are inserted
c     directly into the enlarged valbig array for this piece, if called from bound
c     or entire patch filled if called from filpatch
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

!   ## loop over the 27 regions (in 3D) of the patch - the interior is i=j=k=2 plus
!   ## the ghost cell regions.  If any parts stick out of domain and are periodic
!   ## map indices periodically.  use scratch patch storage of exact size to get rid
!   ## of the iputst, stuff that broke the variable coefficient boundry condition routines.

!   ## if a region sticks out of domain  but is not periodic, not handled in pre-filrecir
!   ## but in bc3amr, now called from bound

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

              mi = i2 - i1 + 1
              mj = j2 - j1 + 1
              mk = k2 - k1 + 1

              if (naux .gt. 0)   
     1            call auxCopyIn(auxPatch,mi,mj,mk,auxbig,
     2                           mitot,mjtot,mktot,naux,
     3                           i1,i2,j1,j2,k1,k2,iglo,jglo,kglo)

              call filrecur(level,nvar,valPatch,auxPatch,naux,time,
     1                      mi,mj,mk,
     2                      1,1,1,                      
     3                      i1+ishift(i), i2+ishift(i),
     4                      j1+jshift(j), j2+jshift(j),
     5                      k1+kshift(k), k2+kshift(k),.true.,msrc)      

                    ! copy it back to proper place in valbig 
                    call patchCopyOut(nvar,valPatch,mi,mj,mk,valbig,
     1                                mitot,mjtot,mktot,
     2                                i1,i2,j1,j2,k1,k2,
     3                                iglo,jglo,kglo)         

           end if

 10      end do
 20     end do
 30    end do

       return
       end
! ============================================================================================

       subroutine patchCopyOut(nvar,valpatch,mi,mj,mk,valbig,
     1                         mitot,mjtot,mktot,
     2                         i1,i2,j1,j2,k1,k2,iglo,jglo,kglo)
 
      ! the patch was filled from a possibly periodically wrapped place.
      ! put it back where it should go in original grids solution array
          
      use amr_module
      implicit none

!     Input
      integer  mi,mj,mk,nvar,mitot,mjtot,mktot
      integer  i1,i2,j1,j2,k1,k2,iglo,jglo,kglo

!     Output
      real*8    valbig(nvar,mitot,mjtot,mktot)
      real*8  valpatch(nvar,mi,mj,mk)

!      Local storage
      integer ist, jst, kst, ivar, i, j, k


      ! this ghost cell patch subset goes from (i1,j1,k1) to (i2,j2,k2) in integer index space
      ! the grid (including ghost cells) is from (iglo,jglo,kglo) to (ighi,jghi,kghi)
      ! figure out where to copy
      ist = i1 - iglo    ! offset by 1 below when coyp, since soln array is 1-based
      jst = j1 - jglo 
      kst = k1 - kglo 

      do k = 1, mk 
      do j = 1, mj 
      do i = 1, mi 
      do ivar = 1, nvar
         valbig(ivar,ist+i,jst+j,kst+k) = valpatch(ivar,i,j,k)
      end do
      end do
      end do
      end do

      return
      end

! ============================================================================================

      subroutine auxCopyIn(auxPatch,mi,mj,mk,auxbig,mitot,mjtot,mktot,
     1                     naux,i1,i2,j1,j2,k1,k2,iglo,jglo,kglo)

      ! set the aux array for the patch  to go with the soln vals to  be filled in filpatch,
      ! by copying from valbig's auxbig array

      use amr_module
      implicit none

!     Input
      integer mi, mj, mk, naux, mitot, mjtot, mktot
      integer i1, i2, j1, j2, k1, k2, iglo, jglo, kglo

!     Output
      real*8  auxbig(naux,mitot,mjtot,mktot)
      real*8  auxPatch(naux,mi,mj,mk)
 
!     Local storage
      integer  ist, jst, kst, iaux , i, j, k


       ! this ghost cell patch subset goes from (i1,j1,k1) to (i2,j2,k2) in integer index space
       ! the grid (including ghost cells) is from (iglo,jglo,kglo) to (ighi,jghi,kghi)
       ! figure out where to copy
       ist = i1 - iglo    ! offset by 1 below  since aux arrays are 1-based
       jst = j1 - jglo 
       kst = k1 - kglo 

      do k = 1, mk 
      do j = 1, mj 
      do i = 1, mi 
      do iaux = 1, naux

       auxPatch(iaux,ist+i,jst+j,kst+k) = auxbig(iaux,i,j,k)
 
      end do
      end do
      end do
      end do

      return
      end
