!> When a fine grid is updated, this subroutine is called to save fluxes
!! going into an adjacent coarse cell in **svdflx**, which is for conservative 
!! fix later.
!!
subroutine fluxad(xfluxm,xfluxp,yfluxm,yfluxp,&
        svdflx,mptr,mitot,mjtot, &
        nvar,lenbc,lratiox,lratioy,ng,dtf,dx,dy)
    !

    use amr_module
    implicit real(CLAW_REAL) (a-h,o-z)


    ! :::::::::::::::::::: FLUXAD ::::::::::::::::::::::::::::::::::
    !  save fine grid fluxes  at the border of the grid, for fixing
    !  up the adjacent coarse cells. at each edge of the grid, only
    !  save the plus or minus fluxes, as necessary. For ex., on
    !  left edge of fine grid, it is the minus xfluxes that modify the
    !  coarse cell.
    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    dimension xfluxm(nvar,mitot,mjtot), yfluxm(nvar,mitot,mjtot)
    dimension xfluxp(nvar,mitot,mjtot), yfluxp(nvar,mitot,mjtot)
    dimension svdflx(nvar,lenbc)

    nx  = mitot-2*ng
    ny  = mjtot-2*ng
    nyc = ny/lratioy
    nxc = nx/lratiox

    ! ::::: left side saved first
    lind = 0

    do j=1,nyc
        lind = lind + 1
        jfine = (j-1)*lratioy + ng
        do ivar = 1, nvar
            do l=1,lratioy
                svdflx(ivar,lind) = svdflx(ivar,lind) + &
                    xfluxm(ivar,ng+1,jfine+l)*dtf*dy
            enddo
        enddo
    enddo

    ! ::::: top side
    !      write(dbugunit,*)" saving top side "
    do i=1,nxc
        lind = lind + 1
        ifine = (i-1)*lratiox + ng
        do ivar = 1, nvar
            do l=1,lratiox
                svdflx(ivar,lind) = svdflx(ivar,lind) + &
                    yfluxp(ivar,ifine+l,mjtot-ng+1)*dtf*dx
            enddo
        enddo
    enddo

    ! ::::: right side
    do j=1,nyc
        lind = lind + 1
        jfine = (j-1)*lratioy + ng
        do ivar = 1, nvar
            do l=1,lratioy
                svdflx(ivar,lind) = svdflx(ivar,lind) + &
                    xfluxp(ivar,mitot-ng+1,jfine+l)*dtf*dy
            enddo
        enddo
    enddo

    ! ::::: bottom side
    !      write(dbugunit,*)" saving bottom side "
    do i=1,nxc
        lind = lind + 1
        ifine = (i-1)*lratiox + ng
        do ivar = 1, nvar
            do l=1,lratiox
                svdflx(ivar,lind) = svdflx(ivar,lind) + &
                    yfluxm(ivar,ifine+l,ng+1)*dtf*dx
            enddo
        enddo
    enddo

    return
end subroutine fluxad
