!  :::::::::::::: PREFILRECUR :::::::::::::::::::::::::::::::::::::::::::
!     For periodic boundary conditions more work needed to fill the
!     piece of the boundary. This routine was
!     called because the patch sticks out of the domain,
!     and has periodic bc.s preprocess the patch before calling
!     filpatch to shift the patch periodically back into the domain.
!
!     Inputs to this routine:
!     xl, xr, yb, yt = the location in physical space of
!     corners of a patch.
!     fill_indices = the location in index space of this patch.
!
!     Outputs from this routine:
!     The values around the border of the grid are inserted
!     directly into the enlarged valbig array for this piece.
!
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
recursive subroutine prefilrecur(level,nvar,valbig,auxbig,naux,time,mitot,  &
                                 nrowst,ilo,ihi,iglo,ighi,patchOnly)



    use amr_module, only: iregsz, nghost, xlower, xperdom
    use amr_module, only: hxposs, NEEDS_TO_BE_SET, alloc
    
    !for setaux timing
    use amr_module, only: timeSetaux, timeSetauxCPU
    
    implicit none

    ! Input
    integer, intent(in) :: level, nvar, naux, mitot
    integer, intent(in) :: ilo,ihi,iglo,ighi
    real(kind=8), intent(in) :: time
    ! false when called from bound, when valbig is whole grid but only filling patch.
    ! true for recursive coarse sub-patches - grid is patch
    logical  :: patchOnly  

    ! Output
    real(kind=8), intent(in out) :: valbig(nvar,mitot)
    real(kind=8), intent(in out) :: auxbig(naux,mitot)
    
    ! Local storage

    ! Various of these are extra. Remove once you know which.
    integer :: i, ii, ivar, ng, i1, i2, nrowst
    integer :: iputst, mi, locpatch, locpaux
    integer :: iwrap1, iwrap2, tmp, locflip, rect(2)
    real(kind=8) :: xlwrap
    integer ::  msrc    ! this signifies not a real grid, no bndry list with it
    ! it is possible to preprocess in the periodic case, just more complicated, so postponing

    integer :: ist(3), iend(3), ishift(3)
    real(kind=8) :: scratch(mitot*nghost*nvar)
    real(kind=8) :: scratchaux(mitot*nghost*naux)

    ! dimension at largest possible
    real(kind=8) :: valPatch((ihi-ilo+1) * nvar)
    real(kind=8) :: auxPatch((ihi-ilo+1) * naux)
    
    !for timing setaux
    integer :: clock_start, clock_finish, clock_rate
    real(kind=8) :: cpu_start, cpu_finish

!     # will divide patch  (from ilo to ihi)  into 3 possibilities (some empty):
!       x sticks out left, x interior, x sticks out right
!       i from (ilo,-1), (0,iregsz(level)-1), (iregsz(level),ihi)
!     # this patch lives in a grid with soln array valbig, which goes from
!       (iglo) to (ighi).

    msrc = -1   ! iitialization indicating whether real src grid so can use faster bc list

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

!   ## loop over the 3 regions (in 1D) of the patch - the interior is i=2 plus
!   ## the ghost cell regions.  If any parts stick out of domain and are periodic
!   ## map indices periodically, but stick the values in the correct place
!   ## in the orig grid (indicated by iputst.
!   ## if a region sticks out of domain  but is not periodic, not handled in (pre)-icall 
!   ## but in setaux/bcamr (not called from here).
    do 20 i = 1, 3
        i1 = max(ilo,  ist(i))
        i2 = min(ihi, iend(i))
        if (i1 .gt. i2) go to 20
            ! make temp patch of just the right size.
            mi = i2 - i1 + 1
            if (mi .gt. (ihi-ilo+1))  then
                write(*,*)" prefilp: not big enough dimension"
            endif
            if (naux .gt. 0)                                                         &
                call auxCopyIn(auxPatch,mi,auxbig,mitot,naux,i1,i2,iglo)

            call filrecur(level,nvar,valPatch,auxPatch,naux,time,mi,1,       &
                    i1+ishift(i),i2+ishift(i),.true.,msrc)

            ! copy it back to proper place in valbig
            call patchCopyOut(nvar,valPatch,mi,valbig,mitot,i1,i2,iglo)

 20 continue

contains

    integer pure function iadd(n,i)
        implicit none
        integer, intent(in) :: n, i
        iadd = locflip + n-1 + nvar*(i-1)
    end function iadd

    integer pure function iaddscratch(n,i)
        implicit none
        integer, intent(in) :: n, i
        iaddscratch = n + nvar*(i-1)  ! no subtract 1

    end function iaddscratch


end subroutine prefilrecur

! ============================================================================================

subroutine patchCopyOut(nvar,valpatch,mi,valbig,mitot,i1,i2,iglo)
 
    ! the patch was filled from a possibly periodically wrapped place.
    ! put it back where it should go in original grids solution array
          
    use amr_module
    implicit none

    ! Input
    integer :: mi, nvar, mitot, i1, i2, iglo, ighi

    ! Output
    real(kind=8), intent(in out) :: valbig(nvar,mitot)
    real(kind=8), intent(in out) :: valpatch(nvar,mi)

    ! Local storage
    integer :: ist


    ! this ghost cell patch subset goes from (i1) to (i2) in integer index space
    ! the grid (including ghost cells) is from (iglo) to (ighi)
    ! figure out where to copy
    ist = i1 - iglo + 1   ! offset 1 since soln array is 1-based

    valbig(:,ist:ist+mi-1) = valpatch

end subroutine patchCopyOut

! ============================================================================================

subroutine auxCopyIn(auxPatch,mi,auxbig,mitot,naux,i1,i2,iglo)

    ! set the aux array for the patch  to go with the soln vals to  be filled in filpatch,
    ! by copying from valbig's auxbig array

    use amr_module
    implicit none

    ! Input
    integer :: mi, naux, mitot, i1, i2, iglo, ighi

    ! Output
    real(kind=8), intent(in out) :: auxbig(naux,mitot)
    real(kind=8), intent(in out) :: auxPatch(naux,mi)

    ! Local storage
    integer :: ist


    ! this ghost cell patch subset goes from (i1,j1) to (i2,j2) in integer index space
    ! the grid (including ghost cells) is from (iglo,jglo) to (ighi,jghi)
    ! figure out where to copy
    ist = i1 - iglo + 1   ! offset 1 since aux arrays are 1-based

    auxPatch(:,1:mi) = auxbig(:,ist:ist+mi-1)

end subroutine auxCopyIn
