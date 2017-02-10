! :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
!
!  fill the portion of valbig from rows  nrowst
!  the patch is described by the corners (xlp) by (xrp).
!  vals are needed at time t, and level level,
!
!  first fill with  values obtainable from the level level
!  grids. if any left unfilled, then enlarge remaining strip of
!  unfilled values by 1 (for later linear interp), and recusively
!  obtain the remaining values from  coarser levels.
!
! :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;
!
recursive subroutine filrecur(level,nvar,valbig,aux,naux,t,mitot, &
     nrowst,ilo,ihi,patchOnly,msrc)

  use amr_module, only: hxposs, xlower, xupper
  use amr_module, only: outunit, nghost, xperdom
  use amr_module, only: iregsz, intratx, NEEDS_TO_BE_SET

  !for setaux timing
  use amr_module, only: timeSetaux, timeSetauxCPU

  implicit none

  ! Input
  integer, intent(in) :: level, nvar, naux, mitot, nrowst
  integer, intent(in) :: ilo,ihi, msrc
  real(kind=8), intent(in) :: t
  logical  :: patchOnly

  ! Output
  real(kind=8), intent(in out) :: valbig(nvar,mitot)
  real(kind=8), intent(in out) :: aux(naux,mitot)

  ! Local storage
  integer :: i_fine, i_coarse, n, k
  integer  :: iplo, iphi
  integer :: mitot_patch, mitot_coarse, nghost_patch, lencrse
  integer :: refinement_ratio_x
  integer :: unset_indices(2)
  real(kind=8) :: dx_fine, dx_coarse
  real(kind=8) :: xlow_coarse, xlow_fine, xhi_fine
  real(kind=8) :: xcent_fine, xcent_coarse,ratiox,floor

  !for timing
  integer :: clock_start, clock_finish, clock_rate
  real(kind=8) :: cpu_start, cpu_finish

  ! Interpolation variables
  real(kind=8) :: eta1, eta2, valp, valm, valc, dupc, dumc
  real(kind=8) :: ducc, du, fu, valint, uslope

  ! Cell set tracking
  logical :: set
  integer(kind=1) :: flaguse(ihi-ilo+1)

  ! Scratch storage
  !  use stack-based scratch arrays instead of alloc, since dont really
  !  need to save beyond these routines, and to allow dynamic memory resizing
  !
  !     use 1d scratch arrays that are potentially the same size as 
  !     current grid, since may not coarsen.
  !     need to make it 1d and do own indexing, since
  !     when pass it in to subroutines they treat it as having dierent
  !     dimensions than the max size need to allocate here

  !--      dimension valcrse((ihi-ilo+2)*nvar)  ! NB this is a 1D array
  !--      dimension auxcrse((ihi-ilo+2)*naux)  ! the +2 is to expand on coarse grid to enclose fine
  ! ### turns out you need 3 rows, forget offset of 1 plus one on each side
  ! the +3 is to expand on coarse grid to enclose fine
  real(kind=8) :: valcrse((ihi-ilo+3)* nvar)
  real(kind=8) :: auxcrse((ihi-ilo+3) * naux)

  ! We begin by filling values for grids at level level. If all values can be
  ! filled in this way, we return;
  mitot_patch = ihi-ilo + 1 ! nrowp

  dx_fine     = hxposs(level)

  ! Coordinates of edges of patch (xlp,xrp)
  xlow_fine = xlower + ilo * dx_fine

  ! Fill in the patch as much as possible using values at this level
  ! note that if only a patch, msrc = -1, otherwise a real grid and intfil
  ! uses its boundary list
  ! msrc either -1 (for a patch) or the real grid number
  call intfil(valbig,mitot,t,flaguse,nrowst,ilo,  &
              ihi,level,nvar,naux,msrc)


  ! Trimbd returns set = true if all of the entries are filled (=1.).
  ! set = false, otherwise. If set = true, then no other levels are
  ! are required to interpolate, and we return.
  !
  ! Note that the used array is filled entirely in intfil, i.e. the
  ! marking done there also takes  into account the points filled by
  ! the boundary conditions. bc2amr will be called later, after all 4
  ! boundary pieces filled.
  call trimbd(flaguse,mitot_patch,set,unset_indices)
  ! il,ir= unset_indices(2)


  ! If set is .true. then all cells have been set and we can skip to setting
  ! the remaining boundary cells.  If it is .false. we need to interpolate
  ! some values from coarser levels, possibly calling this routine
  ! recursively.
  if (.not.set) then

     ! Error check 
     if (level == 1) then
        write(outunit,*)" error in filrecur - level 1 not set"
        write(outunit,'("start at row: ",i4)') nrowst
        print *," error in filrecur - level 1 not set"
        print *," should not need more recursion "
        print *," to set patch boundaries"
        print '("start at row: ",i4)', nrowst
        stop
     endif

     ! We begin by initializing the level level arrays, so that we can use
     ! purely recursive formulation for interpolating.
     dx_coarse  = hxposs(level - 1)

     ! Adjust unset_indices to account for the patch
     ! isl, isr
     unset_indices(1) = unset_indices(1) + ilo - 1
     unset_indices(2) = unset_indices(2) + ilo - 1

     ! Coarsened geometry
     refinement_ratio_x = intratx(level - 1)

     ! New patch rectangle (after we have partially filled it in) but in the
     ! coarse patches [iplo,iphi]

     iplo = (unset_indices(1) - refinement_ratio_x + nghost * refinement_ratio_x) &
          / refinement_ratio_x - nghost
     iphi = (unset_indices(2) + refinement_ratio_x) / refinement_ratio_x

     xlow_coarse = xlower + iplo * dx_coarse

     ! Coarse grid number of spatial points (nrowc)
     mitot_coarse   =  iphi - iplo + 1

     ! Check to make sure we created big enough scratch arrays
     if (mitot_coarse > ihi - ilo + 3) then

        print *," did not make big enough work space in filrecur "
        print *," need coarse space with mitot_coarse,mjtot_coarse ",mitot_coarse
        print *," made space for ilo,ihi ",ilo,ihi
        stop
     endif

     ! Set the aux array values for the coarse grid, this could be done 
     ! instead in intfil using possibly already available bathy data from the
     ! grids
     if (naux > 0) then
        nghost_patch = 0
        lencrse = (ihi-ilo+3)*naux ! set 1 component, not all naux
        do k = 1, lencrse, naux
           auxcrse(k) = NEEDS_TO_BE_SET ! new system checks initialization before setting aux vals
        end do
        call system_clock(clock_start,clock_rate)
        call cpu_time(cpu_start)
        call setaux(nghost_patch, mitot_coarse,  &
             xlow_coarse,            &
             dx_coarse,naux,auxcrse)
        call system_clock(clock_finish,clock_rate)
        call cpu_time(cpu_finish)
        timeSetaux = timeSetaux + clock_finish - clock_start
        timeSetauxCPU = timeSetauxCPU + cpu_finish - cpu_start
     endif

     ! Fill in the edges of the coarse grid. for recursive calls, patch indices and
     ! 'coarse grid' indices are the same (no actual coarse grid here, so cant use mptr
     ! must pass indices. patchOnly argument  is thus true

     if ((xperdom) .and. sticksout(iplo,iphi)) then
        call prefilrecur(level-1,nvar,valcrse,auxcrse,naux,t,mitot_coarse,1,   &
             iplo,iphi,iplo,iphi,.true.)
     else
        call filrecur(level-1,nvar,valcrse,auxcrse,naux,t,mitot_coarse,1,   &
             iplo,iphi,.true.,-1)  ! when going to coarser patch, no source grid (for now at least)
     endif

     ! convert to real for use below
     ratiox = real(refinement_ratio_x,kind=8)

    do i_fine = 1,mitot_patch
        !i_coarse = 2 + (i_fine - (unset_indices(1) - ilo) - 1) / refinement_ratio_x
        !eta1 = (-0.5d0 + real(mod(i_fine - 1, refinement_ratio_x),kind=8)) &
        !                    / real(refinement_ratio_x,kind=8)

        ! new coarse indexing
        !i_coarse =(i_fine+ilo-1)/refinement_ratio_x - iplo + 1

        i_coarse = floor((i_fine+ilo-1)/ratiox) - iplo + 1
        xcent_coarse = xlow_coarse + (i_coarse-.5d0)*dx_coarse
        xcent_fine =  xlower + (i_fine-1+ilo + .5d0)*dx_fine
        eta1 = (xcent_fine-xcent_coarse)/dx_coarse
        if (abs(eta1) .gt. .5) then
            write(*,*)" filpatch x indices wrong: eta1 = ",eta1
        endif

   
        if (flaguse(i_fine) == 0) then

            do n=1,nvar
                ! write(*,*)n,i_coarse+1,coarse_index(n,i_coarse + 1)
                valp = valcrse(coarse_index(n,i_coarse + 1))
                valm = valcrse(coarse_index(n,i_coarse - 1))
                valc   = valcrse(coarse_index(n,i_coarse  ))

                dupc = valp - valc
                dumc = valc   - valm
                ducc = valp - valm
                du   = min(abs(dupc), abs(dumc))
                du   = min(2.d0 * du, 0.5d0 * abs(ducc))
                fu = max(0.d0, sign(1.d0, dupc * dumc))
                uslope = du*sign(1.d0,ducc)*fu ! not really - should divide by h

                valint = valc + eta1 * uslope

                valbig(n,i_fine+nrowst-1) = valint

            end do

        endif
    end do
  end if

  !  set bcs, whether or not recursive calls needed. set any part of patch that stuck out
  xhi_fine = xlower + (ihi + 1) * dx_fine
  ! only call if a small coarser recursive patch
  ! otherwise whole grid bcs done from bound
  if (patchOnly) then
     call bc1amr(valbig,aux,mitot,nvar,naux,dx_fine,level,t,    &
                 xlow_fine,xhi_fine)
  endif

contains

  integer pure function coarse_index(n,i)
    implicit none
    integer, intent(in) :: n, i
    coarse_index = n + nvar*(i-1)
  end function coarse_index

  logical pure function sticksout(iplo,iphi)
    implicit none
    integer, intent(in) :: iplo,iphi
    sticksout = (iplo < 0 .or. iphi >= iregsz(level - 1))
  end function sticksout

end subroutine filrecur
