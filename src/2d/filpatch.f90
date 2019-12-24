! :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
!
!  fill the portion of valbig from rows  nrowst
!                             and  cols  ncolst
!  the patch can also be described by the corners (xlp,ybp) by (xrp,ytp).
!  vals are needed at time t, and level level,
!
!  first fill with  values obtainable from the level level
!  grids. if any left unfilled, then enlarge remaining rectangle of
!  unfilled values by 1 (for later linear interp), and recusively
!  obtain the remaining values from  coarser levels.
!
! :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;
!
!  Below are comments for Doxygen
!> Fill a region (patch) described by:
!! * global indices of its lower left corner: (**ilo**, **jlo**).
!! * global indices of its upper right corner: (**ihi**, **jhi**). 
!!
!! The patch is on grid **msrc** and starts from 
!! row **nrowst** and column **ncolst** of the grid.
!!
!! It first fills the patch with values obtainable from level **level**
!! grids. if any left unfilled, then enlarge remaining rectangle of
!! unfilled values by 1 (for later linear interp), and recusively
!! obtain the remaining values from  coarser levels.
!!
!! ## Algorithm
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!  procedure filrecur(level, patch_fine, val_fine, firstCall)
!!      fill patch_fine as much as possible by copy values from other level "level" grids     
!!      values on patch_fine is stored in array val_fine after the filling
!!      if this is first call to filrecur() (firstCall == true)
!!          don't fill cells outside the computational domain (it will be handled elsewhere)
!!      else
!!          this is a sencond or higher-level (recursive) call to filrecur()
!!          fill cells outside the computational domain by calling bc2amr since they are used for interpolation for finer cells above them
!!      end if
!!      if not all cells in patch_fine is filled (set)
!!          find the smallest rectangular patch on level "level-1", patch_coarse, that contains all unset fine cells in patch_fine
!!          filrecur(level-1, patch_coarse, val_coarse, firstCall=false)
!!          for each unset fine cells in patch_fine
!!              interpolate from val_coarse to fill all unset cells on patch_fine 
!!          end for
!!      end if
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!
!! \param level AMR level the patch is on
!! \param nvar number of equations for the system
!! \param valbig data array for solution \f$q \f$ (cover the whole grid **msrc**)
!! \param aux data array for auxiliary variables 
!! \param naux number of auxiliary variables
!! \param t fill the patch with value at time **t**
!! \param mitot number of cells in *i* direction on grid **msrc** 
!! \param mjtot number of cells in *j* direction on grid **msrc**
!! \param nrowst local *i* index (relative to lower-left corner of grid **msrc**) of the left-most cell of the patch to be filled
!! \param ncolst local *j* index (relative to lower-left corner of grid **msrc**) of the lower-most cell of the patch to be filled
!! \param ilo global *i* index of left-most cell of this patch
!! \param ihi global *i* index of right-most cell of this patch 
!! \param jlo global *j* index of lower-most cell of this patch 
!! \param jhi global *j* index of upper-most cell of this patch 
!! \param patchOnly This is 1) false if this is a first level call to filrecur() and won't set bounday cells outside domain; 2) true
!!if this is a second or higher level (recursive) call to filrecur() and will call bc2amr() to get values for cells outside domain
!!(these cells are only for interpolation to fill cells in first level call to filrecur()
!! \param msrc index of the grid 
!!
!! \callgraph
!! \callergraph


recursive subroutine filrecur(level,nvar,valbig,aux,naux,t,mitot,mjtot, &
     nrowst,ncolst,ilo,ihi,jlo,jhi,patchOnly,msrc,do_aux_copy)

  use amr_module, only: hxposs, hyposs, xlower, ylower, xupper, yupper
  use amr_module, only: outunit, nghost, xperdom, yperdom, spheredom
  use amr_module, only: iregsz, jregsz, intratx, intraty, NEEDS_TO_BE_SET

  implicit none

  ! Input
  integer, intent(in) :: level, nvar, naux, mitot, mjtot, nrowst, ncolst
  integer, intent(in) :: ilo,ihi,jlo,jhi, msrc
  real(kind=8), intent(in) :: t
  logical  :: patchOnly, do_aux_copy,yes_do_aux_copy

  ! Output
  real(kind=8), intent(in out) :: valbig(nvar,mitot,mjtot)
  real(kind=8), intent(in out) :: aux(naux,mitot,mjtot)

  ! Local storage
  integer :: i_fine, j_fine, i_coarse, j_coarse, n, k
  integer  :: iplo, iphi, jplo, jphi
  integer :: mitot_patch, mjtot_patch, mitot_coarse, mjtot_coarse, nghost_patch, lencrse
  integer :: refinement_ratio_x, refinement_ratio_y
  integer :: unset_indices(4)
  real(kind=8) :: dx_fine, dy_fine, dx_coarse, dy_coarse
  real(kind=8) :: xlow_coarse,ylow_coarse, xlow_fine, ylow_fine, xhi_fine,yhi_fine
  real(kind=8) :: xcent_fine, xcent_coarse, ycent_fine, ycent_coarse,ratiox,ratioy,floor    

  !for timing
  integer(kind=8) :: clock_start, clock_finish, clock_rate
  real(kind=8) :: cpu_start, cpu_finish

  ! Interpolation variables
  real(kind=8) :: eta1, eta2, valp10, valm10, valc, valp01, valm01, dupc, dumc
  real(kind=8) :: ducc, du, fu, dvpc, dvmc, dvcc, dv, fv, valint, uslope, vslope

  ! Cell set tracking
  logical :: set
  integer(kind=1) :: flaguse(ihi-ilo+1, jhi-jlo+1)

  ! Scratch storage
  !  use stack-based scratch arrays instead of alloc, since dont really
  !  need to save beyond these routines, and to allow dynamic memory resizing
  !
  !     use 1d scratch arrays that are potentially the same size as 
  !     current grid, since may not coarsen.
  !     need to make it 1d instead of 2 and do own indexing, since
  !     when pass it in to subroutines they treat it as having dierent
  !     dimensions than the max size need to allocate here

  !--      dimension valcrse((ihi-ilo+2)*(jhi-jlo+2)*nvar)  ! NB this is a 1D array 
  !--      dimension auxcrse((ihi-ilo+2)*(jhi-jlo+2)*naux)  ! the +2 is to expand on coarse grid to enclose fine
  ! ### turns out you need 3 rows, forget offset of 1 plus one on each side
  ! the +3 is to expand on coarse grid to enclose fine
  real(kind=8) :: valcrse((ihi-ilo+3) * (jhi-jlo+3) * nvar)  
  real(kind=8) :: auxcrse((ihi-ilo+3) * (jhi-jlo+3) * naux)  

  yes_do_aux_copy = .true.

  ! We begin by filling values for grids at level level. If all values can be
  ! filled in this way, we return;
  mitot_patch = ihi-ilo + 1 ! nrowp
  mjtot_patch = jhi-jlo + 1 

  dx_fine     = hxposs(level)
  dy_fine     = hyposs(level)

  ! Coordinates of edges of patch (xlp,xrp,ybp,ytp)
  xlow_fine = xlower + ilo * dx_fine
  ylow_fine = ylower + jlo * dy_fine

  ! Fill in the patch as much as possible using values at this level
  ! note that if only a patch, msrc = -1, otherwise a real grid and intfil
  ! uses its boundary list
  ! msrc either -1 (for a patch) or the real grid number
  call intfil(valbig,aux,mitot,mjtot,t,flaguse,nrowst,ncolst, ilo,  &
              ihi,jlo,jhi,level,nvar,naux,msrc,do_aux_copy)
 


  ! Trimbd returns set = true if all of the entries are filled (=1.).
  ! set = false, otherwise. If set = true, then no other levels are
  ! are required to interpolate, and we return.
  !
  ! Note that the used array is filled entirely in intfil, i.e. the
  ! marking done there also takes  into account the points filled by
  ! the boundary conditions. bc2amr will be called later, after all 4
  ! boundary pieces filled.
  call trimbd(flaguse,mitot_patch,mjtot_patch,set,unset_indices)
  ! il,ir,jb,jt = unset_indices(4)

  ! If set is .true. then all cells have been set and we can skip to setting
  ! the remaining boundary cells.  If it is .false. we need to interpolate
  ! some values from coarser levels, possibly calling this routine
  ! recursively.
  if (.not.set) then

     ! Error check 
     if (level == 1) then
        write(outunit,*)" error in filrecur - level 1 not set"
        write(outunit,'("start at row: ",i4," col ",i4)') nrowst,ncolst
        print *," error in filrecur - level 1 not set"
        print *," should not need more recursion "
        print *," to set patch boundaries"
        print '("start at row: ",i4," col ",i4)', nrowst,ncolst
        stop
     endif

     ! We begin by initializing the level level arrays, so that we can use
     ! purely recursive formulation for interpolating.
     dx_coarse  = hxposs(level - 1)
     dy_coarse  = hyposs(level - 1)

     ! Adjust unset_indices to account for the patch
     ! isl, isr, jsb, jst
     unset_indices(1) = unset_indices(1) + ilo - 1
     unset_indices(2) = unset_indices(2) + ilo - 1
     unset_indices(3) = unset_indices(3) + jlo - 1
     unset_indices(4) = unset_indices(4) + jlo - 1

     ! Coarsened geometry
     refinement_ratio_x = intratx(level - 1)
     refinement_ratio_y = intraty(level - 1)

     ! New patch rectangle (after we have partially filled it in) but in the
     ! coarse patches [iplo,iphi,jplo,jphi]

     ! islo = unset_indices(1)
     ! ishi = unset_indices(2)
     !                                islo                     ishi 
     !                     |----|----|----|----|----|----|----|----|
     ! If following operation is conducted, it will lead to: 
     !
     !         iplo                                                        iphi
     ! |-------------------|-------------------|-------------------|-------------------|
     iplo = (unset_indices(1) - refinement_ratio_x + nghost * refinement_ratio_x) &
          / refinement_ratio_x - nghost
     iphi = (unset_indices(2) + refinement_ratio_x) / refinement_ratio_x
     jplo = (unset_indices(3) - refinement_ratio_y + nghost * refinement_ratio_y) &
          / refinement_ratio_y - nghost
     jphi = (unset_indices(4) + refinement_ratio_y) / refinement_ratio_y

     xlow_coarse = xlower + iplo * dx_coarse
     ylow_coarse = ylower + jplo * dy_coarse

     ! Coarse grid number of spatial points (nrowc,ncolc)
     mitot_coarse   =  iphi - iplo + 1
     mjtot_coarse   =  jphi - jplo + 1

     ! Check to make sure we created big enough scratch arrays
     if (mitot_coarse > ihi - ilo + 3 .or. &
          mjtot_coarse > jhi - jlo + 3) then

        print *," did not make big enough work space in filrecur "
        print *," need coarse space with mitot_coarse,mjtot_coarse ",mitot_coarse,mjtot_coarse
        print *," made space for ilo,ihi,jlo,jhi ",ilo,ihi,jlo,jhi
        stop
     endif

     ! Set the aux array values for the coarse grid, this could be done 
     ! instead in intfil using possibly already available bathy data from the
     ! grids
     if (naux > 0) then
        nghost_patch = 0
        lencrse = (ihi-ilo+3)*(jhi-jlo+3)*naux ! set 1 component, not all naux
        do k = 1, lencrse, naux
           auxcrse(k) = NEEDS_TO_BE_SET ! new system checks initialization before setting aux vals
        end do
        ! commented out since now copy from coarser grids during intfil
        !call setaux(nghost_patch, mitot_coarse,mjtot_coarse,  &
        !     xlow_coarse, ylow_coarse,            &
        !     dx_coarse,dy_coarse,naux,auxcrse)
     endif

     ! Fill in the edges of the coarse grid. for recursive calls, patch indices and
     ! 'coarse grid' indices are the same (no actual coarse grid here, so cant use mptr
     ! must pass indices. patchOnly argument  is thus true
     if ((xperdom .or. (yperdom .or. spheredom)) .and. sticksout(iplo,iphi,jplo,jphi)) then
        call prefilrecur(level-1,nvar,valcrse,auxcrse,naux,t,mitot_coarse,mjtot_coarse,1,1,   &
             iplo,iphi,jplo,jphi,iplo,iphi,jplo,jphi,.true.)
     else
        call filrecur(level-1,nvar,valcrse,auxcrse,naux,t,mitot_coarse,mjtot_coarse,1,1,   &
             iplo,iphi,jplo,jphi,.true.,-1,yes_do_aux_copy)  ! when going to coarser patch, no source grid (for now at least)
     endif

     ! convert to real for use below
     ratiox = real(refinement_ratio_x,kind=8)
     ratioy = real(refinement_ratio_y,kind=8)

     ! fine below refers to level "level"
     ! coarse below refers to level "level-1"
     do j_fine  = 1,mjtot_patch
        !j_coarse = 2 + (j_fine - (unset_indices(3) - jlo) - 1) / refinement_ratio_y
        !eta2 = (-0.5d0 + real(mod(j_fine - 1, refinement_ratio_y),kind=8)) &
        !                    / real(refinement_ratio_y,kind=8)

        ! we are now processing the j_fine^{th} column of the fine patch
        ! this is contained in the j_coarse^{th} column of the coarse patch below it  
        ! note that these two are not global indices but local indices
        j_coarse =floor((j_fine+jlo-1)/ratioy) - jplo + 1
        ycent_coarse = ylow_coarse + (j_coarse-.5d0)*dy_coarse
        ycent_fine =  ylower + (j_fine-1+jlo + .5d0)*dy_fine
        eta2 = (ycent_fine-ycent_coarse)/dy_coarse
        if (abs(eta2) .gt. .5) then
           write(*,*)" filpatch y indices wrong: eta2 = ",eta2
        endif

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

   
           if (flaguse(i_fine,j_fine) == 0) then

              do n=1,nvar
                 ! write(*,*)n,i_coarse+1,j_coarse,coarse_index(n,i_coarse + 1,j_coarse)
                 ! QUESTION: why do we interpolate in this way?
                 valp10 = valcrse(coarse_index(n,i_coarse + 1,j_coarse))
                 valm10 = valcrse(coarse_index(n,i_coarse - 1,j_coarse))
                 valc   = valcrse(coarse_index(n,i_coarse    ,j_coarse))
                 valp01 = valcrse(coarse_index(n,i_coarse    ,j_coarse + 1))
                 valm01 = valcrse(coarse_index(n,i_coarse    ,j_coarse - 1))

                 dupc = valp10 - valc
                 dumc = valc   - valm10
                 ducc = valp10 - valm10
                 du   = min(abs(dupc), abs(dumc))
                 du   = min(2.d0 * du, 0.5d0 * abs(ducc))
                 fu = max(0.d0, sign(1.d0, dupc * dumc))
                 uslope = du*sign(1.d0,ducc)*fu ! not really - should divide by h

                 dvpc = valp01 - valc
                 dvmc = valc   - valm01
                 dvcc = valp01 - valm01
                 dv   = min(abs(dvpc), abs(dvmc))
                 dv   = min(2.d0 * dv, 0.5d0 * abs(dvcc))
                 fv = max(0.d0,sign(1.d0, dvpc * dvmc))
                 vslope = dv*sign(1.d0,dvcc)*fv ! but faster to put in eta above

                 valint = valc + eta1 * uslope + eta2 * vslope

                 valbig(n,i_fine+nrowst-1,j_fine+ncolst-1) = valint

              end do

           endif
        end do
     end do
  end if  ! end if not set

  !  set bcs, whether or not recursive calls needed. set any part of patch that stuck out
  xhi_fine = xlower + (ihi + 1) * dx_fine
  yhi_fine = ylower + (jhi + 1) * dy_fine
  ! only call if a small coarser recursive patch
  ! otherwise whole grid bcs done from bound
  if (patchOnly) then
     call bc2amr(valbig,aux,mitot,mjtot,nvar,naux,dx_fine,dy_fine,level,t,    &
                 xlow_fine,xhi_fine,ylow_fine,yhi_fine)
  endif

contains

  integer pure function coarse_index(n,i,j)
    implicit none
    integer, intent(in) :: n, i, j
    coarse_index = n + nvar*(i-1)+nvar*mitot_coarse*(j-1)
  end function coarse_index

  logical pure function sticksout(iplo,iphi,jplo,jphi)
    implicit none
    integer, intent(in) :: iplo,iphi,jplo,jphi
    sticksout = (iplo < 0 .or. jplo < 0 .or. &
         iphi >= iregsz(level - 1) .or. jphi >= jregsz(level - 1))
  end function sticksout

end subroutine filrecur
