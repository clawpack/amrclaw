#include "amr_macros.H"

module reflux_module

    use amr_module

    implicit none

    contains

! :::::::::::::::::::: FLUXAD ::::::::::::::::::::::::::::::::::
!  save fine grid fluxes  at the border of the grid, for fixing
!  up the adjacent coarse cells. at each edge of the grid, only
!  save the plus or minus fluxes, as necessary. For ex., on
!  left edge of fine grid, it is the minus xfluxes that modify the
!  coarse cell.
!  We assume this kernel is launched with 2*(mx+my) threads in total,
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! TODO: move this to fluxad.f90
attributes(global) &
subroutine fluxad_gpu(fm, fp, gm, gp, &
    mbc, mx, my, lenbc, lratiox, lratioy, &
    svdflx, dtf, dx, dy)

    real(kind=8), intent(inout) :: fm(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: fp(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: gm(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: gp(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    integer, value, intent(in)  :: mbc, mx, my, lenbc, lratiox, lratioy

    real(kind=8), value, intent(in) :: dtf, dx, dy
    real(kind=8), intent(inout) :: svdflx(NEQNS,lenbc)
    integer :: tid, m
    integer :: i, ifine, j, jfine
    integer :: l
    integer :: mxc, myc

#ifdef DEBUG
    if (blockDim%y /= 1 .or. gridDim%y /= 1) then
        print *, "fluxad_gpu kernel should be called with blockDim%y == 1 and gridDim%y == 1"
        stop
    endif
#endif

    tid = (blockIdx%x-1) * blockDim%x + threadIdx%x

    mxc = mx/2
    myc = my/2

    if (tid <= myc) then
    ! ::::: left side saved first
        j = tid
        jfine = (j-1)*lratioy
        do m = 1,NEQNS
            do l=1,lratioy
                svdflx(m,tid) = svdflx(m,tid) + &
                    fm(1,jfine+l,m)*dtf*dy
            enddo
        enddo

    else if (tid <= myc+mxc) then
    ! ::::: top side
        i = tid - myc
        ifine = (i-1)*lratiox
        do m = 1,NEQNS
            do l=1,lratiox
                svdflx(m,tid) = svdflx(m,tid) + &
                    gp(ifine+l,my+1,m)*dtf*dx
            enddo
        enddo
    else if (tid <= 2*myc+mxc) then
    ! ::::: right side
        j = tid - (myc+mxc)
        jfine = (j-1)*lratioy
        do m = 1,NEQNS
            do l=1,lratioy
                svdflx(m,tid) = svdflx(m,tid) + &
                    fp(mx+1,jfine+l,m)*dtf*dy
            enddo
        enddo
    else if (tid <= 2*(myc+mxc)) then
    ! ::::: bottom side
        i = tid - (2*myc+mxc)
        ifine = (i-1)*lratiox
        do m = 1,NEQNS
            do l=1,lratiox
                svdflx(m,tid) = svdflx(m,tid) + &
                    gm(ifine+l,1,m)*dtf*dx
            enddo
        enddo
    endif
end subroutine fluxad_gpu

attributes(device) &
subroutine fluxad_dev(fm, fp, gm, gp, &
    mbc, mx, my, lenbc, lratiox, lratioy, &
    svdflx, dtf, dx, dy)

    real(kind=8), intent(inout) :: fm(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: fp(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: gm(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    real(kind=8), intent(inout) :: gp(1-mbc:mx+mbc, 1-mbc:my+mbc, NEQNS)
    integer, value, intent(in)  :: mbc, mx, my, lenbc, lratiox, lratioy

    real(kind=8), value, intent(in) :: dtf, dx, dy
    real(kind=8), intent(inout) :: svdflx(NEQNS,lenbc)
    integer :: tid, m
    integer :: i, ifine, j, jfine
    integer :: l
    integer :: mxc, myc

    tid = (blockIdx%x-1) * blockDim%x + threadIdx%x

    mxc = mx/2
    myc = my/2

    if (tid <= myc) then
    ! ::::: left side saved first
        j = tid
        jfine = (j-1)*lratioy
        do m = 1,NEQNS
            do l=1,lratioy
                svdflx(m,tid) = svdflx(m,tid) + &
                    fm(1,jfine+l,m)*dtf*dy
            enddo
        enddo

    else if (tid <= myc+mxc) then
    ! ::::: top side
        i = tid - myc
        ifine = (i-1)*lratiox
        do m = 1,NEQNS
            do l=1,lratiox
                svdflx(m,tid) = svdflx(m,tid) + &
                    gp(ifine+l,my+1,m)*dtf*dx
            enddo
        enddo
    else if (tid <= 2*myc+mxc) then
    ! ::::: right side
        j = tid - (myc+mxc)
        jfine = (j-1)*lratioy
        do m = 1,NEQNS
            do l=1,lratioy
                svdflx(m,tid) = svdflx(m,tid) + &
                    fp(mx+1,jfine+l,m)*dtf*dy
            enddo
        enddo
    else if (tid <= 2*(myc+mxc)) then
    ! ::::: bottom side
        i = tid - (2*myc+mxc)
        ifine = (i-1)*lratiox
        do m = 1,NEQNS
            do l=1,lratiox
                svdflx(m,tid) = svdflx(m,tid) + &
                    gm(ifine+l,1,m)*dtf*dx
            enddo
        enddo
    endif
end subroutine fluxad_dev

attributes(global) &
subroutine fluxad_fused_gpu( &
            grids, fflux,&
            nghost, num_grids, lratiox, lratioy, &
            delt, dx, dy)
    use cuda_module, only: grid_type
    use amr_module
    implicit none

    integer, value, intent(in) :: lratiox, lratioy, num_grids, nghost
    double precision, value, intent(in) :: delt, dx, dy
    type(grid_type), intent(inout) :: grids(num_grids)
    type(gpu_1d_real_ptr_type), intent(in) :: fflux(15000)

    integer :: ng
    integer :: mptr, nx, ny, lenbc

    ng  = (blockIdx%y-1) * blockDim%y + threadIdx%y

    if (ng > num_grids) then
        return
    endif


    mptr = grids(ng)%mptr
    ! No work to do if this is the coarsest level
    if (.not. associated(fflux(mptr)%ptr)) then
        return
    endif
    nx = grids(ng)%nx
    ny = grids(ng)%ny
    lenbc = 2*(nx/lratiox+ny/lratioy)

    ! we call another subroutine to reshape fm ,fp ,gm, gp
    ! and fflux(mptr)%ptr
    call fluxad_dev( &
        grids(ng)%fm, grids(ng)%fp, grids(ng)%gm, grids(ng)%gp, &
        nghost, nx, ny, lenbc, &
        lratiox, lratioy, &
        fflux(mptr)%ptr, delt, dx, dy)
end subroutine fluxad_fused_gpu


attributes(global) &
subroutine fluxsv_gpu(mptr,&
                xfluxm,xfluxp,yfluxm,yfluxp,&
                listbc,&
                fflux,&
                ndimx,ndimy,nvar,maxsp,dtc,hx,hy)

    use amr_module
    implicit none


    integer, value, intent(in) :: mptr
    integer, value, intent(in) :: ndimx, ndimy, nvar, maxsp
    double precision, value, intent(in) :: dtc, hx, hy
    double precision, intent(in) :: xfluxp(ndimx,ndimy,nvar), yfluxp(ndimx,ndimy,nvar)
    double precision, intent(in) :: xfluxm(ndimx,ndimy,nvar), yfluxm(ndimx,ndimy,nvar)
    integer, intent(in) :: listbc(5,maxsp)
    type(gpu_1d_real_ptr_type), intent(in) :: fflux(15000)
    ! local
    integer :: ispot, mkid, intopl, loc
    integer :: i,j, ivar
    integer :: tid

    ! :::::::::::::::::::: FLUXSV :::::::::::::::::::::::::
    !
    !  coarse grids should save their fluxes in cells adjacent to
    !  their nested fine grids, for later conservation fixing.
    !  listbc holds info for where to save which fluxes.
    !  xflux holds 'f' fluxes, yflux holds 'g' fluxes.
    !
    ! :::::::::::::::::::::::::::::;:::::::::::::::::::::::

#ifdef DEBUG
    if (blockDim%y /= 1 .or. gridDim%y /= 1) then
        print *, "fluxsv_gpu kernel should be called with blockDim%y == 1 and gridDim%y == 1"
        stop
    endif
#endif

    tid = (blockIdx%x-1) * blockDim%x + threadIdx%x

    if (tid > maxsp) then 
        return
    endif
    if (listbc(1,tid) .eq. 0) then 
        return
    endif


    ispot = tid

    i        = listbc(1,ispot)
    j        = listbc(2,ispot)
    mkid     = listbc(4,ispot)
    intopl   = listbc(5,ispot)

    loc = nvar*(intopl-1)

    ! side k (listbc 3) has which side of coarse cell has interface
    ! so can save appropriate fluxes.  (dont know why we didnt have
    ! which flux to save directly (i.e. put i+1,j to save that flux
    ! rather than putting in cell center coords).

    if (listbc(3,ispot) .eq. 1) then
        !           ::::: Cell i,j is on right side of a fine grid
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = -xfluxp(i,j,ivar)*dtc*hy
        enddo
    endif

    if (listbc(3,ispot) .eq. 2) then
        !           ::::: Cell i,j on bottom side of fine grid
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = -yfluxm(i,j+1,ivar)*dtc*hx
        enddo
    endif

    if (listbc(3,ispot) .eq. 3) then
        !           ::::: Cell i,j on left side of fine grid
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = -xfluxm(i+1,j,ivar)*dtc*hy
        enddo
    endif

    if (listbc(3,ispot) .eq. 4) then
        !           ::::: Cell i,j on top side of fine grid
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = -yfluxp(i,j,ivar)*dtc*hx 
        enddo
    endif
    !
    !        ### new bcs 5 and 6 come from spherical mapping. note sign change:
    !        ### previous fluxes stored negative flux, fine grids always add
    !        ### their flux, then the delta is either added or subtracted as
    !        ### appropriate for that side.  New bc adds or subtracts BOTH fluxes.
    !
    if (listbc(3,ispot) .eq. 5) then
        !           ::::: Cell i,j on top side of fine grid with spherical mapped bc
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = yfluxm(i,j+1,ivar)*dtc*hx 
        enddo
    endif
    !
    if (listbc(3,ispot) .eq. 6) then
        !           ::::: Cell i,j on bottom side of fine grid with spherical mapped bc
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = yfluxp(i,j,ivar)*dtc*hx
        enddo
    endif
    return
end subroutine fluxsv_gpu

attributes(device) &
subroutine fluxsv_dev(mptr,&
                xfluxm,xfluxp,yfluxm,yfluxp,&
                listbc,&
                fflux,&
                ndimx,ndimy,nvar,maxsp,dtc,hx,hy)

    use amr_module
    implicit none

    integer, value, intent(in) :: mptr
    integer, value, intent(in) :: ndimx, ndimy, nvar, maxsp
    double precision, value, intent(in) :: dtc, hx, hy
    double precision, intent(in) :: xfluxp(ndimx,ndimy,nvar), yfluxp(ndimx,ndimy,nvar)
    double precision, intent(in) :: xfluxm(ndimx,ndimy,nvar), yfluxm(ndimx,ndimy,nvar)
    integer, intent(in) :: listbc(5,maxsp)
    type(gpu_1d_real_ptr_type), intent(in) :: fflux(15000)
    ! local
    integer :: ispot, mkid, intopl, loc
    integer :: i,j, ivar
    integer :: tid

    ! :::::::::::::::::::: FLUXSV :::::::::::::::::::::::::
    !
    !  coarse grids should save their fluxes in cells adjacent to
    !  their nested fine grids, for later conservation fixing.
    !  listbc holds info for where to save which fluxes.
    !  xflux holds 'f' fluxes, yflux holds 'g' fluxes.
    !
    ! :::::::::::::::::::::::::::::;:::::::::::::::::::::::

    tid = (blockIdx%x-1) * blockDim%x + threadIdx%x

    if (tid > maxsp) then 
        return
    endif
    if (listbc(1,tid) .eq. 0) then 
        return
    endif

    ispot = tid

    i        = listbc(1,ispot)
    j        = listbc(2,ispot)
    mkid     = listbc(4,ispot)
    intopl   = listbc(5,ispot)

    loc = nvar*(intopl-1)

    ! side k (listbc 3) has which side of coarse cell has interface
    ! so can save appropriate fluxes.  (dont know why we didnt have
    ! which flux to save directly (i.e. put i+1,j to save that flux
    ! rather than putting in cell center coords).

    if (listbc(3,ispot) .eq. 1) then
        !           ::::: Cell i,j is on right side of a fine grid
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = -xfluxp(i,j,ivar)*dtc*hy
        enddo
    endif

    if (listbc(3,ispot) .eq. 2) then
        !           ::::: Cell i,j on bottom side of fine grid
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = -yfluxm(i,j+1,ivar)*dtc*hx
        enddo
    endif

    if (listbc(3,ispot) .eq. 3) then
        !           ::::: Cell i,j on left side of fine grid
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = -xfluxm(i+1,j,ivar)*dtc*hy
        enddo
    endif

    if (listbc(3,ispot) .eq. 4) then
        !           ::::: Cell i,j on top side of fine grid
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = -yfluxp(i,j,ivar)*dtc*hx 
        enddo
    endif
    !
    !        ### new bcs 5 and 6 come from spherical mapping. note sign change:
    !        ### previous fluxes stored negative flux, fine grids always add
    !        ### their flux, then the delta is either added or subtracted as
    !        ### appropriate for that side.  New bc adds or subtracts BOTH fluxes.
    !
    if (listbc(3,ispot) .eq. 5) then
        !           ::::: Cell i,j on top side of fine grid with spherical mapped bc
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = yfluxm(i,j+1,ivar)*dtc*hx 
        enddo
    endif
    !
    if (listbc(3,ispot) .eq. 6) then
        !           ::::: Cell i,j on bottom side of fine grid with spherical mapped bc
        do ivar = 1, nvar
            fflux(mkid)%ptr(loc + ivar) = yfluxp(i,j,ivar)*dtc*hx
        enddo
    endif
    return
end subroutine fluxsv_dev

attributes(global) &
subroutine fluxsv_fused_gpu(grids, cflux, fflux, &
    nghost,num_grids,nvar,maxsp,delt,hx,hy)

    use amr_module
    use cuda_module, only: grid_type
    implicit none

    integer, value, intent(in) :: nghost, nvar, maxsp, num_grids
    double precision, value, intent(in) :: delt,hx,hy
    type(grid_type), intent(inout) :: grids(num_grids)
    type(gpu_1d_real_ptr_type), intent(in) :: fflux(15000)
    type(gpu_2d_int_ptr_type), intent(in)  :: cflux(15000)

    integer :: mptr, ng, nx, ny, mitot, mjtot

    ng  = (blockIdx%y-1) * blockDim%y + threadIdx%y

    if (ng > num_grids) then
        return
    endif

    mptr = grids(ng)%mptr
    ! No work to do if this is the finest level
    if (.not. associated(cflux(mptr)%ptr)) then
        return
    endif
    nx = grids(ng)%nx
    ny = grids(ng)%ny
    mitot  = nx + 2*nghost
    mjtot  = ny + 2*nghost

    ! we call another subroutine to reshape fm ,fp ,gm, gp etc.
    call fluxsv_dev(mptr, &
             grids(ng)%fm,grids(ng)%fp, grids(ng)%gm, grids(ng)%gp, &
             cflux(mptr)%ptr, &
             fflux, &
             mitot,mjtot,nvar,maxsp,delt,hx,hy)
end subroutine fluxsv_fused_gpu


!> For each coarse-fine interface, a Riemann problem between an inner
!! ghost cell value on the fine grid and cell value in the adjacent coarse
!! cell must be solved and added to corresponding location in
!! **node(ffluxptr, mptr)** for conservative fix later
!!
! -------------------------------------------------------------
!
subroutine qad_cpu(valbig,mitot,mjtot,nvar, &
        svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy,&
        maux,aux,auxc1d,delt,mptr)

    use amr_module
#ifdef PROFILE
    use profiling_module
#endif


    !
    ! ::::::::::::::::::::::::::: QAD ::::::::::::::::::::::::::::::::::
    !  are added in to coarse grid value, as a conservation fixup. 
    !  Done each fine grid time step. If source terms are present, the
    !  coarse grid value is advanced by source terms each fine time step too.

    !  No change needed in this sub. for spherical mapping: correctly
    !  mapped vals already in bcs on this fine grid and coarse saved
    !  vals also properly prepared
    !
    ! Side 1 is the left side of the fine grid patch.  Then go around clockwise.
    ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    !
    !      # local storage
    !      # note that dimension here are bigger than dimensions used
    !      # in rp2, but shouldn't matter since wave is not used in qad
    !      # and for other arrays it is only the last parameter that is wrong
    !      #  ok as long as meqn, mwaves < maxvar

    integer, parameter :: max1dp1 = max1d+1
    integer :: mitot, mjtot, nvar, lenbc, maux
    integer :: lratiox, lratioy, mptr
    integer :: iaddaux, iaux, nc, nr, level, index
    integer :: ivar, lind, ncrse
    double precision :: hx, hy, delt, tgrid
    double precision :: ql(nvar,max1dp1),    qr(nvar,max1dp1)
    double precision :: wave(nvar,mwaves,max1dp1), s(mwaves,max1dp1)
    double precision :: amdq(nvar,max1dp1),  apdq(nvar,max1dp1)
    double precision :: auxl(maxaux*max1dp1),  auxr(maxaux*max1dp1)
    double precision :: valbig(nvar,mitot,mjtot)
    double precision :: qc1d(nvar,lenbc)
    double precision :: svdflx(nvar,lenbc)
    double precision :: aux(maux,mitot,mjtot)
    double precision :: auxc1d(maux,lenbc)

    logical :: qprint

    integer :: i,j,ma,ic,jc,l,influx,ifine,jfine

    !
    !  WARNING: auxl,auxr dimensioned at max possible, but used as if
    !  they were dimensioned as the real maux by max1dp1. Would be better
    !  of course to dimension by maux by max1dp1 but this wont work if maux=0
    !  So need to access using your own indexing into auxl,auxr.
    iaddaux(iaux,i) = iaux + maux*(i-1)

    data qprint/.false./
    !
    !      aux is auxiliary array with user parameters needed in Riemann solvers
    !          on fine grid corresponding to valbig
    !      auxc1d is coarse grid stuff from around boundary, same format as qc1d
    !      auxl, auxr are work arrays needed to pass stuff to rpn2
    !      maux is the number of aux variables, which may be zero.
    !

#ifdef PROFILE
    call nvtxStartRange("qad",13)
#endif
    tgrid = rnode(timemult, mptr)
    if (qprint) &
        write(dbugunit,*)" working on grid ",mptr," time ",tgrid
    nc = mjtot-2*nghost
    nr = mitot-2*nghost
    level = node(nestlevel, mptr)
    index = 0

    !
    !--------
    !  side 1
    !--------
    !
    do j = nghost+1, mjtot-nghost
        if (maux.gt.0) then
            do ma = 1,maux
                if (auxtype(ma).eq."xleft") then
                    !                # Assuming velocity at left-face, this fix
                    !                # preserves conservation in incompressible flow:
                    auxl(iaddaux(ma,j-nghost+1)) = aux(ma,nghost+1,j)
                else
                    !                # Normal case -- we set the aux arrays 
                    !                # from the cell corresponding  to q
                    auxl(iaddaux(ma,j-nghost+1)) = aux(ma,nghost,j)
                endif
            enddo
        endif
        do ivar = 1, nvar
            ql(ivar,j-nghost+1) = valbig(ivar,nghost,j)
        enddo
    enddo

    lind = 0
    ncrse = (mjtot-2*nghost)/lratioy
    do jc = 1, ncrse
        index = index + 1
        do l = 1, lratioy
            lind = lind + 1
            if (maux.gt.0) then
                do ma=1,maux
                    auxr(iaddaux(ma,lind)) = auxc1d(ma,index)
                enddo
            endif
            do ivar = 1, nvar
                qr(ivar,lind) = qc1d(ivar,index)
            enddo
        enddo
    enddo

    if (qprint) then
        write(dbugunit,*) 'side 1, ql and qr:'
        do i=2,nc
            write(dbugunit,4101) i,qr(1,i-1),ql(1,i)
        enddo
        4101      format(i3,4e16.6)
        if (maux .gt. 0) then
            write(dbugunit,*) 'side 1, auxr:'
            do i=2,nc
                write(dbugunit,4101) i,(auxr(iaddaux(ma,i-1)),ma=1,maux)
            enddo
            write(dbugunit,*) 'side 1, auxl:'
            do i=2,nc
                write(dbugunit,4101) i,(auxl(iaddaux(ma,i)),ma=1,maux)
            enddo
        endif
    endif

    call rpn2(1,max1dp1-2*nghost,nvar,mwaves,maux,nghost, &
        nc+1-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)
    !
    ! we have the wave. for side 1 add into sdflxm
    !
    influx = 0
    do j = 1, nc/lratioy
        influx  = influx + 1
        jfine = (j-1)*lratioy
        do ivar = 1, nvar
            do l = 1, lratioy
                svdflx(ivar,influx) = svdflx(ivar,influx) &
                    + amdq(ivar,jfine+l+1) * hy * delt &
                    + apdq(ivar,jfine+l+1) * hy * delt
            enddo
        enddo
    enddo

    !--------
    !  side 2
    !--------
    !
    if (mjtot .eq. 2*nghost+1) then
        !          # a single row of interior cells only happens when using the
        !          # 2d amrclaw code to do a 1d problem with refinement.
        !          # (feature added in Version 4.3)
        !          # skip over sides 2 and 4 in this case
        go to 299
    endif

    do i = nghost+1, mitot-nghost
        if (maux.gt.0) then
            do ma = 1,maux
                auxr(iaddaux(ma,i-nghost)) = aux(ma,i,mjtot-nghost+1)
            enddo
        endif
        do ivar = 1, nvar
            qr(ivar,i-nghost) = valbig(ivar,i,mjtot-nghost+1)
        enddo
    enddo

    lind = 0
    ncrse = (mitot-2*nghost)/lratiox
    do ic = 1, ncrse
        index = index + 1
        do l = 1, lratiox
            lind = lind + 1
            if (maux.gt.0) then
                do ma=1,maux
                    if (auxtype(ma).eq."yleft") then
                        !                # Assuming velocity at bottom-face, this fix
                        !                # preserves conservation in incompressible flow:
                        ifine = (ic-1)*lratiox + nghost + l
                        auxl(iaddaux(ma,lind+1)) = aux(ma,ifine,mjtot-nghost+1)
                    else
                        auxl(iaddaux(ma,lind+1)) = auxc1d(ma,index)
                    endif
                enddo
            endif
            do ivar = 1, nvar
                ql(ivar,lind+1) = qc1d(ivar,index)
            enddo
        enddo
    enddo

    if (qprint) then
        write(dbugunit,*) 'side 2, ql and qr:'
        do i=1,nr
            write(dbugunit,4101) i,ql(1,i+1),qr(1,i)
        enddo
        if (maux .gt. 0) then
            write(dbugunit,*) 'side 2, auxr:'
            do i = 1, nr
                write(dbugunit,4101) i, (auxr(iaddaux(ma,i)),ma=1,maux)
            enddo
            write(dbugunit,*) 'side 2, auxl:'
            do i = 1, nr
                write(dbugunit,4101) i, (auxl(iaddaux(ma,i)),ma=1,maux)
            enddo
        endif
    endif
    call rpn2(2,max1dp1-2*nghost,nvar,mwaves,maux,nghost, &
        nr+1-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)
    !
    ! we have the wave. for side 2. add into sdflxp
    !
    do i = 1, nr/lratiox
        influx  = influx + 1
        ifine = (i-1)*lratiox
        do ivar = 1, nvar
            do l = 1, lratiox
                svdflx(ivar,influx) = svdflx(ivar,influx) &
                    - amdq(ivar,ifine+l+1) * hx * delt &
                    - apdq(ivar,ifine+l+1) * hx * delt
            enddo
        enddo
    enddo

    299  continue

    !--------
    !  side 3
    !--------
    !
    do j = nghost+1, mjtot-nghost
        if (maux.gt.0) then
            do ma = 1,maux
                auxr(iaddaux(ma,j-nghost)) = aux(ma,mitot-nghost+1,j)
            enddo
        endif
        do ivar = 1, nvar
            qr(ivar,j-nghost) = valbig(ivar,mitot-nghost+1,j)
        enddo
    enddo

    lind = 0
    ncrse = (mjtot-2*nghost)/lratioy
    do jc = 1, ncrse
        index = index + 1
        do l = 1, lratioy
            lind = lind + 1
            if (maux.gt.0) then
                do ma=1,maux
                    if (auxtype(ma).eq."xleft") then
                        !                # Assuming velocity at left-face, this fix
                        !                # preserves conservation in incompressible flow:
                        jfine = (jc-1)*lratioy + nghost + l
                        auxl(iaddaux(ma,lind+1)) = aux(ma,mitot-nghost+1,jfine)
                    else
                        auxl(iaddaux(ma,lind+1)) = auxc1d(ma,index)
                    endif
                enddo
            endif
            do ivar = 1, nvar
                ql(ivar,lind+1) = qc1d(ivar,index)
            enddo
        enddo
    enddo

    if (qprint) then
        write(dbugunit,*) 'side 3, ql and qr:'
        do i=1,nc
            write(dbugunit,4101) i,ql(1,i),qr(1,i)
        enddo
    endif
    call rpn2(1,max1dp1-2*nghost,nvar,mwaves,maux,nghost, &
        nc+1-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)
    !
    ! we have the wave. for side 3 add into sdflxp
    !
    do j = 1, nc/lratioy
        influx  = influx + 1
        jfine = (j-1)*lratioy
        do ivar = 1, nvar
            do l = 1, lratioy
                svdflx(ivar,influx) = svdflx(ivar,influx) &
                    - amdq(ivar,jfine+l+1) * hy * delt &
                    - apdq(ivar,jfine+l+1) * hy * delt
            enddo
        enddo
    enddo

    !--------
    !  side 4
    !--------
    !
    if (mjtot .eq. 2*nghost+1) then
        !          # a single row of interior cells only happens when using the
        !          # 2d amrclaw code to do a 1d problem with refinement.
        !          # (feature added in Version 4.3)
        !          # skip over sides 2 and 4 in this case
        go to 499
    endif
    !
    do i = nghost+1, mitot-nghost
        if (maux.gt.0) then
            do ma = 1,maux
                if (auxtype(ma).eq."yleft") then
                    !                # Assuming velocity at bottom-face, this fix
                    !                # preserves conservation in incompressible flow:
                    auxl(iaddaux(ma,i-nghost+1)) = aux(ma,i,nghost+1)
                else
                    auxl(iaddaux(ma,i-nghost+1)) = aux(ma,i,nghost)
                endif
            enddo
        endif
        do ivar = 1, nvar
            ql(ivar,i-nghost+1) = valbig(ivar,i,nghost)
        enddo
    enddo

    lind = 0
    ncrse = (mitot-2*nghost)/lratiox
    do ic = 1, ncrse
        index = index + 1
        do l = 1, lratiox
            lind = lind + 1
            if (maux.gt.0) then
                do ma=1,maux
                    auxr(iaddaux(ma,lind)) = auxc1d(ma,index)
                enddo
            endif
            do ivar = 1, nvar
                qr(ivar,lind) = qc1d(ivar,index)
            enddo
        enddo
    enddo

    if (qprint) then
        write(dbugunit,*) 'side 4, ql and qr:'
        do i=1,nr
            write(dbugunit,4101) i, ql(1,i),qr(1,i)
        enddo
    endif
    call rpn2(2,max1dp1-2*nghost,nvar,mwaves,maux,nghost, &
        nr+1-2*nghost,ql,qr,auxl,auxr,wave,s,amdq,apdq)
    !
    ! we have the wave. for side 4. add into sdflxm
    !
    do i = 1, nr/lratiox
        influx  = influx + 1
        ifine = (i-1)*lratiox
        do ivar = 1, nvar
            do l = 1, lratiox
                svdflx(ivar,influx) = svdflx(ivar,influx) &
                    + amdq(ivar,ifine+l+1) * hx * delt &
                    + apdq(ivar,ifine+l+1) * hx * delt
            enddo
        enddo
    enddo

    499   continue

    !      # for source terms:
    if (method(5) .ne. 0) then   ! should I test here if index=0 and all skipped?
        call src1d(nvar,nghost,lenbc,qc1d,maux,auxc1d,tgrid,delt)
        !      # how can this be right - where is the integrated src term used?
    endif

#ifdef PROFILE
    call nvtxEndRange()
#endif
    return
end subroutine qad_cpu

end module reflux_module
