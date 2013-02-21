! :::::::::::::::::::::::::::::: FILVAL ::::::::::::::::::::::::::
!
! create and fill coarser (lev-1) patch with one extra coarse cell all
! around, plus the ghost cells . will interpolate from this patch to grid mptr 
! without needing special boundary code. 
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine filval(val,mitot,mjtot,hx,hy,lev,time,valc,auxc,mic,mjc,xleft, &
                  xright,ybot,ytop,nvar,mptr,ilo,ihi,jlo,jhi,aux,naux,locflip)
 
    use amr_module, only: intratx, intraty, max1d, xperdom, yperdom, spheredom
    use amr_module, only: xlower, xupper, ylower, yupper, mcapa, nghost
    implicit none

    ! Input
    integer, intent(in) :: mitot, mjtot, lev, mic, mjc, nvar, mptr, ilo, ihi
    integer, intent(in) :: jlo, jhi, naux, locflip
    real(kind=8), intent(in) :: hx, hy, time, xleft, xright, ybot, ytop

    ! Output
    real(kind=8), intent(in out) :: val(nvar,mitot,mjtot), valc(nvar,mic,mjc)
    real(kind=8), intent(in out) :: aux(naux,mitot,mjtot), auxc(naux,mic,mjc)

    ! Locals
    integer :: levc, lratiox, lratioy, iclo, jclo, ichi, jchi, ng, i, j, ivar
    integer :: ico, jco, ifine, jfine
    real(kind=8) :: hxcrse, hycrse, xl, xr, yb, yt, slp, slm, slopex, slopey
    real(kind=8) :: xoff, yoff
    
    real(kind=8) :: dudx(max1d), dudy(max1d)

    levc    = lev - 1
    lratiox = intratx(levc)
    lratioy = intraty(levc)
    hxcrse  = hx*lratiox
    hycrse  = hy*lratioy
    xl      = xleft  - hxcrse 
    xr      = xright + hxcrse
    yb      = ybot   - hycrse 
    yt      = ytop   + hycrse

    ! set integer indices for coarser patch enlarged by 1 cell 
    ! (can stick out of domain). proper nesting will insure this one
    ! call is sufficient.
    iclo   = ilo/lratiox - 1
    jclo   = jlo/lratioy - 1
    ichi   = (ihi+1)/lratiox - 1 + 1
    jchi   = (jhi+1)/lratioy - 1 + 1
    ng     = 0

    ! :::  mcapa  is the capacity function index
    if (mcapa == 0) then   !dont need to copy aux stuff along with soln
        if (xperdom .or. yperdom .or. spheredom) then
            call preintcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,locflip)
        else
            call intcopy(valc,mic,mjc,nvar,iclo,ichi,jclo,jchi,levc,1,1)
        endif
    else  ! intersect grids and copy all (soln and aux)
        if (xperdom .or. yperdom .or. spheredom) then
            call preicall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi, &
                          levc,locflip)
        else
            call icall(valc,auxc,mic,mjc,nvar,naux,iclo,ichi,jclo,jchi,levc,1,1)
        endif
    endif
    call bc2amr(valc,auxc,mic,mjc,nvar,naux,hxcrse,hycrse,levc,time,xl,xr,yb, &
                yt,xlower,ylower,xupper,yupper,xperdom,yperdom,spheredom)

    ! prepare slopes - use min-mod limiters

    do j = 2, mjc-1
        do ivar = 1, nvar
            do i = 2, mic-1

                slp = valc(ivar,i+1,j) - valc(ivar,i,j)
                slm = valc(ivar,i,j)   - valc(ivar,i-1,j)
                slopex = min(abs(slp), abs(slm))  &
                            * sign(1.0d0, valc(ivar,i+1,j) - valc(ivar,i-1,j))

                ! if there's a sign change, set slope to 0.
                if ( slm*slp .gt. 0.d0) then
                    dudx(i) = slopex
                else
                    dudx(i) = 0.d0
                endif

                slp = valc(ivar,i,j+1) - valc(ivar,i,j)
                slm = valc(ivar,i,j)   - valc(ivar,i,j-1)
                slopey = min(abs(slp), abs(slm))  &
                            * sign(1.0d0, valc(ivar,i,j+1) - valc(ivar,i,j-1))
                if ( slm*slp .gt. 0.d0) then
                    dudy(i) = slopey
                else
                    dudy(i) = 0.d0
                endif

            enddo

            ! interp. from coarse cells to fine grid
            do ico = 1,lratiox
                xoff = (float(ico) - .5)/lratiox - .5
                do jco = 1,lratioy
                    jfine = (j-2)*lratioy + nghost + jco
                    yoff  = (float(jco) - .5)/lratioy - .5
                    do i = 2, mic-1
                        ifine   = (i-2)*lratiox + nghost + ico
                        val(ivar,ifine,jfine) = valc(ivar,i,j) + xoff*dudx(i) &
                                                               + yoff*dudy(i)
                    enddo
                enddo
            enddo
        enddo
    enddo

    if (mcapa .ne. 0) then
        call fixcapaq(val,aux,mitot,mjtot,valc,auxc,mic,mjc,nvar,naux,levc)
    endif

    ! overwrite interpolated values with fine grid values, if available.
    call intcopy(val,mitot,mjtot,nvar,ilo-nghost,ihi+nghost,jlo-nghost, &
                 jhi+nghost,lev,1,1)

end subroutine filval
