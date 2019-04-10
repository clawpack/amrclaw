!> Compute enclosing rectangle around flagged points.
!! save some info., even tho. usage might be low and rect. scrapped.
!! 
!! \param[in] badpts  x,y coords of cells around which a rectangle 
!! is to be made and compute efficiency on
!! \param[in] npt     num. of badpts. in the cluster.
!! \param[out] usage   ratio of flagged to unflagged badpts. in new grid
!!                measures goodness of fit and clustering
!! \param[out] intrect stores some info. for grid created herein.
!!                sometimes rect = rnode, sometimes = temp. array.
!!                sometimes intrect = node.
!!                depending on calling prog. (grdfit or expand)

subroutine moment(intrect, badpts, npt, usage)
    !
    ! :::::::::::::::::::::::: MOMENT ::::::::::::::::::::::::::::::::::
    !  moment = compute enclosing rectangle around flagged points.
    !  save some info., even tho. usage might be low and rect. scrapped.
    !
    ! input parameters:
    !     badpts      = x,y coords of flagged badpts grouped into clusters
    !                   are in the first two rows
    !     npt         = num. of badpts. in the cluster.
    !
    ! output parameters:
    !     usage       = ratio of flagged to unflagged badpts. in new grid
    !                   measures goodness of fit and clustering
    !    intrect( )    = stores some info. for grid created herein.
    !                   sometimes rect = rnode, sometimes = temp. array.
    !                   sometimes intrect = node.
    !                   depending on calling prog. (grdfit or expand)
    !
    !
    ! :::::::::::::::::::::::: MOMENT ::::::::::::::::::::::::::::::::::
    !

    use amr_module, only: ndilo, ndjlo, ndihi, ndjhi, nsize

    ! Input
    integer, intent(in) :: npt
    real(kind=8), intent(in) :: badpts(2, npt)
    integer, intent(in out) ::  intrect(nsize)
    real(kind=8), intent(out) :: usage

    ! Locals
    real(kind=8) :: e_min(2), e_max(2), side(2)

    ! Compute length of enclosing rectangles to include all flagged badpts
    do i = 1, 2
        e_min(i) = minval(badpts(i, :))
        e_max(i) = maxval(badpts(i, :))
    end do

    ! From length of the sides determine rect corners
    ! transform to cell numebrs (subtract 0.5)
    intrect(ndilo) = nint(e_min(1) - 0.5d0)
    intrect(ndjlo) = nint(e_min(2) - 0.5d0)
    intrect(ndihi) = nint(e_max(1) - 0.5d0)
    intrect(ndjhi) = nint(e_max(2) - 0.5d0)

    ! Compute usage
    side(1) = real(intrect(ndihi) - intrect(ndilo) + 1, kind=8)
    side(2) = real(intrect(ndjhi) - intrect(ndjlo) + 1, kind=8)
    usage = real(npt, kind=8) / (side(1) * side(2))

end subroutine moment
