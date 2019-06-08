!> Compute enclosing line around flagged points.
!! save some info., even tho. usage might be low and rect. scrapped.
!! 
!! \param[in] badpts  x coords of cells around which a rectangle 
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

    use amr_module, only: ndilo, ndihi, nsize

    ! Input
    integer, intent(in) :: npt
    real(kind=8), intent(in) :: badpts(1, npt)
    integer, intent(in out) ::  intrect(nsize)
    real(kind=8), intent(out) :: usage

    ! Locals
    real(kind=8) :: e_min, e_max, side

    ! Compute length of enclosing lines to include all flagged badpts
    e_min = minval(badpts(1, :))
    e_max = maxval(badpts(1, :))

    ! From length of the sides determine rect corners
    ! transform to cell numebrs (subtract 0.5)
    intrect(ndilo) = nint(e_min - 0.5d0)
    intrect(ndihi) = nint(e_max - 0.5d0)

    ! Compute usage
    side = real(intrect(ndihi) - intrect(ndilo) + 1, kind=8)
    usage = real(npt, kind=8) / side

end subroutine moment
