!  :::::::::::::: PREFILRECUR :::::::::::::::::::::::::::::::::::::::::::
!>     For periodic boundary conditions more work needed to fill the
!!     piece of the boundary. This routine was
!!     called because the patch sticks out of the domain,
!!     and has periodic bc.s preprocess the patch before calling
!!     filpatch to shift the patch periodically back into the domain.
!!
!!     Inputs to this routine:
!!     xl, xr, yb, yt = the location in physical space of
!!     corners of a patch.
!!     fill_indices = the location in index space of this patch.
!!
!!     Outputs from this routine:
!!     The values around the border of the grid are inserted
!!     directly into the enlarged valbig array for this piece.
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
recursive subroutine prefilrecur(level,nvar,valbig,auxbig,naux,time,mitot,mjtot,    &  
                                 nrowst,ncolst,ilo,ihi,jlo,jhi,iglo,ighi,jglo,jghi, &
                                 patchOnly)



    use amr_module, only: iregsz, jregsz, nghost, xlower, ylower, xperdom, yperdom
    use amr_module, only: spheredom, hxposs, hyposs, NEEDS_TO_BE_SET, alloc
    
    implicit none

    ! Input
    integer, intent(in) :: level, nvar, naux, mitot, mjtot
    integer, intent(in) :: ilo,ihi,jlo,jhi,iglo,ighi,jglo,jghi
    real(kind=8), intent(in) :: time
    ! false when called from bound, when valbig is whole grid but only filling patch.
    ! true for recursive coarse sub-patches - grid is patch
    logical  :: patchOnly  
    logical :: do_aux_copy

    ! Output
    real(kind=8), intent(in out) :: valbig(nvar,mitot,mjtot)
    real(kind=8), intent(in out) :: auxbig(naux,mitot,mjtot)
    
    ! Local storage
    integer :: i, j, ii, jj, ivar, ng, i1, i2, j1, j2, nrowst, ncolst
    integer :: iuse1, iuse2
    integer :: iputst, jputst, mi, mj, locpatch, locpaux
    integer :: jbump, iwrap1, iwrap2, jwrap1, tmp, locflip, rect(4)
    real(kind=8) :: xlwrap, ybwrap
    integer ::  msrc    ! this signifies not a real grid, no bndry list with it
    ! it is possible to preprocess in the periodic case, just more complicated, so postponing

    integer :: ist(3), iend(3), jst(3), jend(3), ishift(3), jshift(3)
    real(kind=8) :: scratch(max(mitot,mjtot)*nghost*nvar)
    real(kind=8) :: scratchaux(max(mitot,mjtot)*nghost*naux)

    ! dimension at largest possible
    !real(kind=8) :: valPatch((ihi-ilo+1) * (jhi-jlo+1) * nvar)  
    !real(kind=8) :: auxPatch((ihi-ilo+1) * (jhi-jlo+1) * naux)  
    real(kind=8) :: valPatch((ihi-ilo+2) * (jhi-jlo+2) * nvar)  
    real(kind=8) :: auxPatch((ihi-ilo+2) * (jhi-jlo+2) * naux)  
    

!     # will divide patch  (from ilo,jlo to ihi,jhi)  into 9 possibilities (some empty): 
!       x sticks out left, x interior, x sticks out right
!       same for y. for example, the max. would be
!       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)
!     # this patch lives in a grid with soln array valbig, which goes from
!       (iglo,jglo) to (ighi,jghi).

    msrc = -1   ! iitialization indicating whether real src grid so can use faster bc list
    do_aux_copy = .false.  ! too complicated in periodic case, just call setaux

    if (xperdom) then       
       ist(1)    = ilo
       ist(2)    = 0
       ist(3)    = iregsz(level)
       iend(1)   = -1
       iend(2)   = iregsz(level)-1
       iend(3)   = ihi
       ishift(1) = iregsz(level)
       ishift(2) = 0
       ishift(3) = -iregsz(level)
    else  ! if not periodic, set vals to only have nonnull intersection for interior regoin
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

    if (yperdom .or. spheredom) then
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
       !jst(1)    = jregsz(level)
       !jst(2)    = jlo
       !jst(3)    = jregsz(level)
       !jend(1)   = -1
       !jend(2)   = jhi
       !jend(3)   = -1
       !jshift(1) = 0
       !jshift(2) = 0
       !jshift(3) = 0
    endif

!   ## loop over the 9 regions (in 2D) of the patch - the interior is i=j=2 plus
!   ## the ghost cell regions.  If any parts stick out of domain and are periodic
!   ## map indices periodically, but stick the values in the correct place
!   ## in the orig grid (indicated by iputst,jputst.
!   ## if a region sticks out of domain  but is not periodic, not handled in (pre)-icall 
!   ## but in setaux/bcamr (not called from here).
!!    if (.not. (xperdom .and. yperdom)) go to 66
    do 20 i = 1, 3
        i1 = max(ilo,  ist(i))
        i2 = min(ihi, iend(i))
        if (i1 .gt. i2) go to 20

    do 10 j = 1, 3
       if (yperdom) then
         j1 = max(jlo,  jst(j))
         j2 = min(jhi, jend(j))
       else ! not periodic in y
          ! initialize to skip then reset below
          j1 = jregsz(level) + nghost + 1
          j2 = -nghost-1
          jshift = 0
          if (j .eq. 1) then
             if (jlo .le. -1) then 
                 j1 = jlo  
                 j2 = 0 ! (top ghost cell + 1 for extrap bc)
             endif
          !endif
          elseif (j .eq. 2) then 
             ! interior part of patch,no need to add 1
             j1 = max(jlo,0)
             j2 = min(jhi,jregsz(level)-1)
          !endif
          elseif (j .eq. 3) then 
             if (jhi .gt. jregsz(level)) then
               j1 = jregsz(level)-1 ! includes last cell on top for extrap bc
               j2 = jhi
             endif
           endif
       endif

        ! part of patch in this region
         if (j1 <= j2) then 
             if (.not. spheredom) then
             !if (.not. spheredom .or. j .eq. 2) then
                 ! make temp patch of just the right size. 
                 mi = i2 - i1 + 1
                 mj = j2 - j1 + 1
                 !if (mi .gt. (ihi-ilo+1) .or.  mj .gt. (jhi-jlo+1))  then
                 if (mi .gt. (ihi-ilo+2) .or.  mj .gt. (jhi-jlo+2))  then
                    write(*,*)" prefilp: not big enough dimension"
                 endif
                 if (naux .gt. 0)                                                         &
                     call auxCopyIn(auxPatch,mi,mj,auxbig,mitot,mjtot,naux,i1,i2,j1,j2,   &
                                    iglo,jglo)
                                       
                 call filrecur(level,nvar,valPatch,auxPatch,naux,time,mi,mj,       &
                         1,1,i1+ishift(i),i2+ishift(i),j1+jshift(j),j2+jshift(j),  &
                               .true.,msrc,do_aux_copy)
                  ! copy it back to proper place in valbig 
                   call patchCopyOut(nvar,valPatch,mi,mj,valbig,mitot,mjtot,i1,i2,j1,j2,   &
                                     iglo,jglo)
 
             else
                mi = i2 - i1 + 1
                mj = j2 - j1 + 1
                ng = 0    ! no ghost cells in this little patch. fill everything.

                jbump = 0
                if (j1 < 0)   jbump = abs(j1)  ! bump up so new bottom is 0
                if (j2 >= jregsz(level)) jbump = -(j2+1-jregsz(level)) ! bump down

                ! next 2 lines would take care of periodicity in x
                iwrap1 = i1 + ishift(i)
                iwrap2 = i2 + ishift(i)
                ! next 2 lines take care of mapped sphere bcs
                iwrap1 = iregsz(level) - iwrap1 -1
                iwrap2 = iregsz(level) - iwrap2 -1
                ! swap so that smaller one is left index, etc since mapping reflects
                tmp = iwrap1
                iwrap1 = iwrap2
                iwrap2 = tmp

                jwrap1 = j1 + jbump
                xlwrap = xlower + iwrap1*hxposs(level)
                ybwrap = ylower + jwrap1*hyposs(level)

                if (naux>0) then
                    scratchaux = NEEDS_TO_BE_SET  !flag all cells with signal since dimensioned strangely
                    call setaux(ng,mi,mj,xlwrap,ybwrap,hxposs(level),hyposs(level),naux,scratchaux)
                endif 

                rect = [iwrap1,iwrap2,j1+jbump,j2+jbump]
                call filrecur(level,nvar,scratch,scratchaux,naux,time,mi, &
                              mj,1,1,iwrap1,iwrap2,j1+jbump,j2+jbump,.false.,msrc,  &
                              do_aux_copy)

                ! copy back using weird mapping for spherical folding (so cant call copy subr below)
                do ii = i1, i2
                do jj = j1, j2
                       ! write(dbugunit,'(" filling loc ",2i5," with ",2i5)') & 
                       ! nrowst+ii-fill_indices(1),ncolst+jj-fill_indices(3),mi-(ii-i1),mj-jj+j1

                       do ivar = 1, nvar
                           valbig(ivar,nrowst+(ii-ilo),ncolst+(jj-jlo)) = &
                               scratch(iaddscratch(ivar,mi-(ii-i1),mj-(jj-j1)))
                       end do
                       ! write(dbugunit,'(" new val is ",4e15.7)')(valbig(ivar,  &
                       ! nrowst+(ii-fill_indices(1)),ncolst+(jj-fill_indices(3))),ivar=1,nvar)
                end do
                end do 
             endif ! end if not spherical or j == 2
         endif ! end if region not empty

  10     continue
 20 continue
 
 return
           
contains

    integer pure function iadd(n,i,j)
        implicit none
        integer, intent(in) :: n, i, j
        iadd = locflip + n-1 + nvar*((j-1)*mi+i-1)
    end function iadd

    integer pure function iaddscratch(n,i,j)
        implicit none
        integer, intent(in) :: n, i, j
        iaddscratch = n + nvar*((j-1)*mi+i-1)  ! no subtract 1

    end function iaddscratch


end subroutine prefilrecur

! ============================================================================================

subroutine patchCopyOut(nvar,valpatch,mi,mj,valbig,mitot,mjtot,i1,i2,j1,j2,iglo,jglo)
 
    ! the patch was filled from a possibly periodically wrapped place.
    ! put it back where it should go in original grids solution array
          
    use amr_module
    implicit none

    ! Input
    integer :: mi, mj, nvar, mitot, mjtot, i1, i2,j1, j2, iglo, ighi, jglo, jghi

    ! Output
    real(kind=8), intent(in out) :: valbig(nvar,mitot,mjtot)
    real(kind=8), intent(in out) :: valpatch(nvar,mi,mj)

    ! Local storage
    integer :: ist, jst 


    ! this ghost cell patch subset goes from (i1,j1) to (i2,j2) in integer index space
    ! the grid (including ghost cells) is from (iglo,jglo) to (ighi,jghi)
    ! figure out where to copy
    ist = i1 - iglo + 1   ! offset 1 since soln array is 1-based
    jst = j1 - jglo + 1

    valbig(:,ist:ist+mi-1, jst:jst+mj-1) = valpatch

end subroutine patchCopyOut

! ============================================================================================

subroutine auxCopyIn(auxPatch,mi,mj,auxbig,mitot,mjtot,naux,i1,i2,j1,j2,iglo,jglo)

    ! set the aux array for the patch  to go with the soln vals to  be filled in filpatch,
    ! by copying from valbig's auxbig array

    use amr_module
    implicit none

    ! Input
    integer :: mi, mj, naux, mitot, mjtot, i1, i2,j1, j2, iglo, ighi, jglo, jghi

    ! Output
    real(kind=8), intent(in out) :: auxbig(naux,mitot,mjtot)
    real(kind=8), intent(in out) :: auxPatch(naux,mi,mj)

    ! Local storage
    integer :: ist, jst 


    ! this ghost cell patch subset goes from (i1,j1) to (i2,j2) in integer index space
    ! the grid (including ghost cells) is from (iglo,jglo) to (ighi,jghi)
    ! figure out where to copy
    ist = i1 - iglo + 1   ! offset 1 since aux arrays are 1-based
    jst = j1 - jglo + 1

    auxPatch(:,1:mi,1:mj) = auxbig(:,ist:ist+mi-1, jst:jst+mj-1)

end subroutine auxCopyIn
