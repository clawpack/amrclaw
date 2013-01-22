!
!  :::::::::::::: BOUND :::::::::::::::::::::::::::::::::::::::::::
!     This routine sets the boundary values for a given grid 
!     at level level.
!     We are setting the values for a strip ng zones wide all
!     the way around the border, in 4 rectangular strips.
!
!     Outputs from this routine:
!     The values around the border of the grid are inserted
!     directly into the enlarged valbig array.
!
!     This routine calls the routine filpatch
!     which for any block of mesh points on a given level,
!     intersects that block with all grids on that level and with
!     the physical boundaries, copies the values into the
!     appropriate intersecting regions, and interpolates the remaining
!     cells from coarser grids as required.
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine bound(time,nvar,ng,valbig,mitot,mjtot,mptr,aux,naux)

    use amr_module, only: rnode, node, hxposs, hyposs, cornxlo, cornxhi
    use amr_module, only: cornylo, cornyhi, ndilo, ndihi, ndjlo, ndjhi
    use amr_module, only: nestlevel, xlower, xupper, ylower, yupper
    use amr_module, only: xperdom, yperdom, spheredom
    
    implicit none

    ! Input
    integer, intent(in) :: nvar, ng, mitot, mjtot, mptr, naux
    real(kind=8), intent(in) :: time
    real(kind=8), intent(in out) :: valbig(nvar,mitot,mjtot)
    real(kind=8), intent(in out) :: aux(naux,mitot,mjtot)

    ! Locals
    integer :: ilo, ihi, jlo, jhi, level, rect(4),i,j
    real(kind=8) :: xleft, xright, ybot, ytop, hx, hy, xl, xr, yb, yt

    xleft  = rnode(cornxlo, mptr)
    xright = rnode(cornxhi, mptr)
    ybot   = rnode(cornylo, mptr)
    ytop   = rnode(cornyhi, mptr)
    ilo    = node(ndilo, mptr)
    ihi    = node(ndihi, mptr)
    jlo    = node(ndjlo, mptr)
    jhi    = node(ndjhi, mptr)
    level  = node(nestlevel, mptr)
    hx     = hxposs(level)
    hy     = hyposs(level)


    ! left boundary
    xl = xleft - ng*hx
    xr = xleft
    yb = ybot 
    yt = ytop

    rect = [ilo-ng,ilo-1,jlo,jhi]
    if ((xl < xlower) .and. xperdom) then
        call  prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,1,ng+1, &
                          rect)
    else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,1,ng+1,rect)
    endif

    ! right boundary
    xl = xright
    xr = xright + ng*hx
    yb = ybot
    yt = ytop

    rect = [ihi+1,ihi+ng,jlo,jhi]
    if ((xr .gt. xupper) .and. xperdom) then
        call  prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
                          mitot-ng+1,ng+1,rect)
    else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
                      mitot-ng+1,ng+1,rect)
    endif


    ! bottom boundary
    xl = xleft  - ng*hx
    xr = xright + ng*hx
    yb = ybot - ng*hy
    yt = ybot
    
    rect = [ilo-ng,ihi+ng,jlo-ng,jlo-1]
    if ( ((yb < ylower) .and. (yperdom .or. spheredom)) .or. &
        ( ((xl < xlower) .or. (xr > xupper)) .and. xperdom) ) then
       call prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,1,1,rect)
    else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,1,1,rect)
    endif

    ! top boundary
    xl = xleft - ng*hx
    xr = xright + ng*hx
    yb = ytop
    yt = ytop + ng*hy

    rect = [ilo-ng,ihi+ng,jhi+1,jhi+ng]
    if ( ((yt .gt. yupper) .and. (yperdom .or. spheredom)) .or. &
         (((xl .lt. xlower) .or. (xr .gt. xupper)) .and. xperdom) ) then
        call prefilrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
                         1,mjtot-ng+1,rect)
    else
        call filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot, &
                      1,mjtot-ng+1,rect)
    endif

! external boundary conditions   THIS IS NOW DONE IN THE FILPATCHES 
! (in the recursive filrecur.f, since filpatches had to call bc2amr, have them 
! all do it).

end subroutine bound
