! :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
!
!  fill the portion of valbig from rows  nrowst
!                             and  cols  ncolst
!  the patch can also be described by the corners (xlp,ybp) by (xrp,ytp).
!  vals are needed at time time , and level level,
!
!  first fill with  values obtainable from the level level
!  grids. if any left unfilled, then enlarge remaining rectangle of
!  unfilled values by 1 (for later linear interp), and recusively
!  obtain the remaining values from  coarser levels.
!
! :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;
recursive subroutine filrecur(level,nvar,valbig,aux,naux,time,mitot,mjtot,nrowst,ncolst,ilo,ihi,jlo,jhi)

    use amr_module, only: hxposs, hyposs, xlower, ylower, xupper, yupper
    use amr_module, only: outunit, nghost, xperdom, yperdom, spheredom
    use amr_module, only: iregsz, jregsz, intratx, intraty

    implicit none

    ! Input
    integer, intent(in) :: level, nvar, naux, mitot, mjtot, nrowst, ncolst
    integer, intent(in) :: ilo, ihi, jlo, jhi
    real(kind=8), intent(in) :: time

    ! Output
    real(kind=8), intent(in out) :: valbig(nvar,mitot,mjtot)
    real(kind=8), intent(in out) :: aux(naux,mitot,mjtot)

    ! Local storage
    integer :: nrowp, ncolp, levc, isl, isr, jsb, jst, lratiox, lratioy, iplo
    integer :: jplo, iphi, jphi, nrowc, ncolc, ntot, iff, jf, ic, jc, ivar
    integer :: il, ir, jb, jt
    real(kind=8) :: hxf, hyf, xlp, xrp, ybp, ytp, hxc, hyc, xlc, ybc, xrc, ytc
    real(kind=8) :: eta1, eta2, valp10, valm10, valc, valp01, valm01, dupc, dumc
    real(kind=8) :: ducc, du, fu, dvpc, dvmc, dvcc, dv, fv, valint
    logical :: set

    ! Scratch storage
    !  use stack-based scratch arrays instead of alloc, since dont really
    !  need to save beyond these routines, and to allow dynamic memory resizing
    !
    !     use 1d scratch arrays that are potentially the same size as 
    !     current grid, since may not coarsen.
    !     need to make it 1d instead of 2 and do own indexing, since
    !     when pass it in to subroutines they treat it as having different
    !     dimensions than the max size need to allocate here
    
    !--      dimension valcrse((ihi-ilo+2)*(jhi-jlo+2)*nvar)  ! NB this is a 1D array 
    !--      dimension auxcrse((ihi-ilo+2)*(jhi-jlo+2)*naux)  ! the +2 is to expand on coarse grid to enclose fine
    ! ### turns out you need 3 rows, forget offset of 1 plus one on each side
    real(kind=8) :: valcrse((ihi-ilo+3)*(jhi-jlo+3)*nvar)  ! NB this is a 1D array 
    real(kind=8) :: auxcrse((ihi-ilo+3)*(jhi-jlo+3)*naux)  ! the +3 is to expand on coarse grid to enclose fine
    integer(kind=1) :: flaguse(ihi-ilo+1,jhi-jlo+1)


    ! We begin by filling values for grids at level level. If all values can be
    ! filled in this way, we return;
    nrowp   = ihi - ilo + 1
    ncolp   = jhi - jlo + 1

    hxf     = hxposs(level)
    hyf     = hyposs(level)
    xlp     = xlower + ilo*hxf
    xrp     = xlower + (ihi+1)*hxf
    ybp     = ylower + jlo*hyf
    ytp     = ylower + (jhi+1)*hyf

    call intfil(valbig,mitot,mjtot,time,flaguse,nrowst,ncolst,ilo,ihi,jlo,jhi,level,nvar,naux)

    ! Trimbd returns set = true if all of the entries are filled (=1.).
    ! set = false, otherwise. If set = true, then no other levels are
    ! are required to interpolate, and we return.
    !
    ! Note that the used array is filled entirely in intfil, i.e. the
    ! marking done there also takes  into account the points filled by
    ! the boundary conditions. bc2amr will be called later, after all 4
    ! boundary pieces filled.
    call trimbd(flaguse,nrowp,ncolp,set,il,ir,jb,jt)

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
        levc = level - 1
        hxc  = hxposs(levc)
        hyc  = hyposs(levc)

        isl  = il + ilo - 1
        isr  = ir + ilo - 1
        jsb  = jb + jlo - 1
        jst  = jt + jlo - 1

        ! Coarsened geometry
        lratiox = intratx(levc)
        lratioy = intraty(levc)
        iplo   = (isl-lratiox+nghost*lratiox)/lratiox - nghost
        jplo   = (jsb-lratioy+nghost*lratioy)/lratioy - nghost
        iphi   = (isr+lratiox)/lratiox
        jphi   = (jst+lratioy)/lratioy

        xlc  =  xlower + iplo*hxc
        ybc  =  ylower + jplo*hyc
        xrc  =  xlower + (iphi+1)*hxc
        ytc  =  ylower + (jphi+1)*hyc

        nrowc   =  iphi - iplo + 1
        ncolc   =  jphi - jplo + 1
        ntot    = nrowc*ncolc*(nvar+naux)
!        write(*,'(" needed coarse grid size ",2i5," allocated ",2i5)') nrowc,ncolc, ihi-ilo+2,jhi-jlo+2
!        write(*,'(" needed coarse grid size ",2i5," allocated ",2i5)') nrowc,ncolc, ihi-ilo+3,jhi-jlo+3
        if (nrowc .gt. ihi-ilo+3 .or. ncolc .gt. jhi-jlo+3) then
            print *," did not make big enough work space in filrecur "
            print *," need coarse space with nrowc,ncolc ",nrowc,ncolc
            print *," made space for ilo,ihi,jlo,jhi ",ilo,ihi,jlo,jhi
            stop
        endif

        ! Set the aux array values, this could be done instead in intfil using
        ! possibly already available bathy data from the grids
        if (naux > 0) then
            call setaux(nrowc - 2*nghost,ncolc - 2*nghost,nghost, &
                        nrowc - 2*nghost,ncolc - 2*nghost, &
                        xlc + nghost*hxc,ybc + nghost*hyc, &
                        hxc,hyc,naux,auxcrse)
        endif

        ! *** NEED TO LOOK INTO THIS CALL AND WHY IT IS DONE ***
        if ((xperdom .or. (yperdom .or. spheredom)) .and. sticksout(iplo,iphi,jplo,jphi)) then
            call prefilrecur(levc,nvar,valcrse,auxcrse,naux,time,nrowc,ncolc,1,1,iplo,iphi,jplo,jphi)
        else
            call filrecur(levc,nvar,valcrse,auxcrse,naux,time,nrowc,ncolc,1,1,iplo,iphi,jplo,jphi)
        endif

        do iff = 1,nrowp
            ic = 2 + (iff - (isl - ilo) - 1) / lratiox
            eta1 = (-0.5d0 + real(mod(iff-1, lratiox),kind=8)) &
                                / real(lratiox,kind=8)

            do jf  = 1,ncolp
                jc = 2 + (jf - (jsb - jlo) - 1) / lratioy
                eta2 = (-0.5d0 + real(mod(jf -1, lratioy),kind=8)) &
                                    / real(lratioy,kind=8)

                if (flaguse(iff,jf) == 0) then

                    do ivar = 1,nvar

                        valp10 = valcrse(ivalc(ivar,ic+1,jc))
                        valm10 = valcrse(ivalc(ivar,ic-1,jc))
                        valc   = valcrse(ivalc(ivar,ic  ,jc))
                        valp01 = valcrse(ivalc(ivar,ic  ,jc+1))
                        valm01 = valcrse(ivalc(ivar,ic  ,jc-1))
        
                        dupc = valp10 - valc
                        dumc = valc   - valm10
                        ducc = valp10 - valm10
                        du   = dmin1(dabs(dupc),dabs(dumc))
                        du   = dmin1(2.d0*du,.5d0*dabs(ducc))
                        fu = dmax1(0.d0,dsign(1.d0,dupc*dumc))
        
                        dvpc = valp01 - valc
                        dvmc = valc   - valm01
                        dvcc = valp01 - valm01
                        dv   = dmin1(dabs(dvpc),dabs(dvmc))
                        dv   = dmin1(2.d0*dv,.5d0*dabs(dvcc))
                        fv = dmax1(0.d0,dsign(1.d0,dvpc*dvmc))

                        valint = valc + eta1 * du * sign(1.d0, ducc) * fu &
                                      + eta2 * dv * sign(1.d0, dvcc) * fv


                        valbig(ivar,iff+nrowst-1,jf+ncolst-1) = valint

                    end do
                endif
            end do
        end do
    end if

    !  set bcs, whether or not recursive calls needed. set any part of patch that stuck out
    call bc2amr(valbig,aux,mitot,mjtot,nvar,naux,hxf,hyf,level,time,xlp,xrp, &
                ybp,ytp,xlower,ylower,xupper,yupper,xperdom,yperdom,spheredom)

contains

    integer pure function ivalc(ivar,i,j)
        implicit none
        integer, intent(in) :: ivar, i, j
        ivalc = ivar + nvar*(i-1)+nvar*nrowc*(j-1)
    end function ivalc

    logical pure function sticksout(iplo,iphi,jplo,jphi)
        implicit none
        integer, intent(in) :: iplo, iphi, jplo, jphi
        sticksout = (iplo < 0 .or. jplo < 0 .or. &
                     iphi >= iregsz(levc) .or. jphi >= jregsz(levc))
    end function sticksout

end subroutine filrecur