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



! Note that valbig here is in SoA format
subroutine qad_cpu2(valbig,mitot,mjtot,mbc,nvar, &
        svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy,&
        delt,mptr,cc,zz)

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

    implicit none

    integer, value, intent(in) :: mitot, mjtot, mbc, nvar, lenbc
    integer, value, intent(in) :: lratiox, lratioy, mptr
    double precision, value, intent(in) :: hx, hy, delt, cc, zz
    double precision, intent(in) :: valbig(mitot,mjtot,nvar)
    double precision, intent(in) :: qc1d(nvar,lenbc)
    double precision, intent(inout) :: svdflx(nvar,lenbc)

    double precision :: ql(nvar,max1d,SPACEDIM*2), qr(nvar,max1d,SPACEDIM*2)
    double precision :: wave(nvar,mwaves,max1d,SPACEDIM*2), s(mwaves,max1d,SPACEDIM*2)
    double precision :: amdq(nvar,max1d,SPACEDIM*2),  apdq(nvar,max1d,SPACEDIM*2)

    integer :: nc, nr
    integer :: ivar, ncrse
    integer :: i,j,ic,jc,l
    integer :: base
    integer :: tid
    double precision :: add

#ifdef PROFILE
    call nvtxStartRange("qad",13)
#endif
    nc = mjtot-2*mbc
    nr = mitot-2*mbc
    do tid = 1, (nc+nr)*2

        ! prepare ql and qr
        !
        ! side 1
        base = 0
        if (tid > base .and. tid <= (base+nc)) then
            j = tid - base
            do ivar = 1, nvar
                ql(ivar,j,1) = valbig(mbc,j+mbc,ivar)
            enddo

            jc = (j-1)/lratioy+1
            do ivar = 1, nvar
                qr(ivar,j,1) = qc1d(ivar,jc)
            enddo
        endif

        ! side 2
        base = nc
        if (tid > base .and. tid <= (base+nr)) then
            i = tid - base
            do ivar = 1, nvar
                qr(ivar,i,2) = valbig(i+mbc,mjtot-mbc+1,ivar)
            enddo

            ic = (i-1)/lratiox+1 + nc/lratioy
            do ivar = 1, nvar
                ql(ivar,i,2) = qc1d(ivar,ic)
            enddo
        endif

        ! side 3
        base = nc + nr
        if (tid > base .and. tid <= (base+nc)) then
            j = tid - base
            do ivar = 1, nvar
                qr(ivar,j,3) = valbig(mitot-mbc+1,j+mbc,ivar)
            enddo

            jc = (j-1)/lratioy+1 + nc/lratioy + nr/lratiox
            do ivar = 1, nvar
                ql(ivar,j,3) = qc1d(ivar,jc)
            enddo
        endif

        ! side 4
        base = 2*nc + nr
        if (tid > base .and. tid <= (base+nr)) then
            i = tid - base
            do ivar = 1, nvar
                ql(ivar,i,4) = valbig(i+mbc,mbc,ivar)
            enddo

            ic = (i-1)/lratiox+1 + 2*nc/lratioy + nr/lratiox
            do ivar = 1, nvar
                qr(ivar,i,4) = qc1d(ivar,ic)
            enddo
        endif

        ! We only need amdq and apdq from this
        call rpn2_all_edges(max1d,mbc, &
                ql,qr,amdq,apdq,tid,nc,nr,cc,zz)

    ! Need to make sure all tasks above are finished before
    ! Starting below
    ! But we can't synchronize arocss CUDA blocks (before CUDA 9.0)
    ! So we must rely on the execution order of CUDA threads
        !--------
        !  side 1
        !--------
        ! we have the wave. for side 1 add into sdflxm
        !
        base = 0
        if (tid > base .and. tid <= (base+nc)) then
            j = tid-base
            jc = (j-1)/lratioy + 1
            do ivar = 1, nvar
                add = amdq(ivar,j,1) * hy * delt &
                    + apdq(ivar,j,1) * hy * delt
                svdflx(ivar,jc) = svdflx(ivar,jc) + add
            enddo
        endif

        !--------
        !  side 2
        !--------
        ! we have the wave. for side 2. add into sdflxp
        !
        base = nc
        if (tid > base .and. tid <= (base+nr)) then
            i = tid-base
            ic = (i-1)/lratiox + 1 + nc/lratioy
            do ivar = 1, nvar
                add = - amdq(ivar,i,2) * hx * delt &
                    - apdq(ivar,i,2) * hx * delt
                svdflx(ivar,ic) = svdflx(ivar,ic) + add
            enddo
        endif

        !--------
        !  side 3
        !--------
        ! we have the wave. for side 3 add into sdflxp
        !
        base = nc + nr
        if (tid > base .and. tid <= (base+nc)) then
            j = tid-base
            jc = (j-1)/lratioy + 1 + nc/lratioy + nr/lratiox
            do ivar = 1, nvar
                add = - amdq(ivar,j,3) * hy * delt &
                    - apdq(ivar,j,3) * hy * delt
                svdflx(ivar,jc) = svdflx(ivar,jc) + add
            enddo
        endif

        !--------
        !  side 4
        !--------
        ! we have the wave. for side 4. add into sdflxm
        !
        base = 2*nc + nr
        if (tid > base .and. tid <= (base+nr)) then
            i = tid-base
            ic = (i-1)/lratiox + 1 + 2*nc/lratioy + nr/lratiox
            do ivar = 1, nvar
                add = amdq(ivar,i,4) * hx * delt &
                    + apdq(ivar,i,4) * hx * delt
                svdflx(ivar,ic) = svdflx(ivar,ic) + add
            enddo
        endif
    enddo

#ifdef PROFILE
    call nvtxEndRange()
#endif
    return
end subroutine qad_cpu2

subroutine rpn2_all_edges(maxm,mbc,ql,qr,amdq,apdq,tid,nc,nr,cc,zz)

    use amr_module

    implicit none

    integer, value, intent(in) :: maxm, mbc, tid, nc, nr
    double precision, value, intent(in) :: cc,zz
    double precision, intent(in)  ::   ql(NEQNS, 1:maxm,SPACEDIM*2)
    double precision, intent(in)  ::   qr(NEQNS, 1:maxm,SPACEDIM*2)
    double precision, intent(out) :: apdq(NEQNS, 1:maxm,SPACEDIM*2)
    double precision, intent(out) :: amdq(NEQNS, 1:maxm,SPACEDIM*2)

    double precision :: wave(NEQNS, NWAVES)
    double precision ::    s(1:NWAVES)

!     local arrays
!     ------------
    double precision :: delta1, delta2, delta3, a1, a2
    integer :: mu, mv, i, m, iedge, base

    ! side 1
    base = 0
    if (tid > base .and. tid <= (base+nc)) then
        iedge = 1
        mu = 2
        mv = 3
        i = tid - base
        delta1 = ql( 1,i,iedge) - qr( 1,i,iedge)
        delta2 = ql(mu,i,iedge) - qr(mu,i,iedge)
        a1 = (-delta1 + zz*delta2) / (2.d0*zz)
        a2 = (delta1 + zz*delta2) / (2.d0*zz)
    
    !        # Compute the waves.
    
        wave( 1,1) = -a1*zz
        wave(mu,1) = a1
        wave(mv,1) = 0.d0
        s(1) = -cc
    
        wave( 1,2) = a2*zz
        wave(mu,2) = a2
        wave(mv,2) = 0.d0
        s(2) = cc
        do m = 1,NEQNS
            amdq(m,i,iedge) = s(1)*wave(m,1)
            apdq(m,i,iedge) = s(2)*wave(m,2)
        enddo
    endif

    ! side 2
    base = nc
    if (tid > base .and. tid <= (base+nr)) then
        iedge = 2
        mu = 3
        mv = 2
        i = tid - base
        delta1 = ql( 1,i,iedge) - qr( 1,i,iedge)
        delta2 = ql(mu,i,iedge) - qr(mu,i,iedge)
        a1 = (-delta1 + zz*delta2) / (2.d0*zz)
        a2 = (delta1 + zz*delta2) / (2.d0*zz)
    
    !        # Compute the waves.
    
        wave( 1,1) = -a1*zz
        wave(mu,1) = a1
        wave(mv,1) = 0.d0
        s(1) = -cc
    
        wave( 1,2) = a2*zz
        wave(mu,2) = a2
        wave(mv,2) = 0.d0
        s(2) = cc
        do m = 1,NEQNS
            amdq(m,i,iedge) = s(1)*wave(m,1)
            apdq(m,i,iedge) = s(2)*wave(m,2)
        enddo
    endif

    ! side 3
    base = nc + nr
    if (tid > base .and. tid <= (base+nc)) then
        iedge = 3
        mu = 2
        mv = 3
        i = tid - base
        delta1 = ql( 1,i,iedge) - qr( 1,i,iedge)
        delta2 = ql(mu,i,iedge) - qr(mu,i,iedge)
        a1 = (-delta1 + zz*delta2) / (2.d0*zz)
        a2 = (delta1 + zz*delta2) / (2.d0*zz)
    
    !        # Compute the waves.
    
        wave( 1,1) = -a1*zz
        wave(mu,1) = a1
        wave(mv,1) = 0.d0
        s(1) = -cc
    
        wave( 1,2) = a2*zz
        wave(mu,2) = a2
        wave(mv,2) = 0.d0
        s(2) = cc
        do m = 1,NEQNS
            amdq(m,i,iedge) = s(1)*wave(m,1)
            apdq(m,i,iedge) = s(2)*wave(m,2)
        enddo
    endif

    ! side 4
    base = 2*nc + nr
    if (tid > base .and. tid <= (base+nr)) then
        iedge = 4
        mu = 3
        mv = 2
        i = tid - base
        delta1 = ql( 1,i,iedge) - qr( 1,i,iedge)
        delta2 = ql(mu,i,iedge) - qr(mu,i,iedge)
        a1 = (-delta1 + zz*delta2) / (2.d0*zz)
        a2 = (delta1 + zz*delta2) / (2.d0*zz)
    
    !        # Compute the waves.
    
        wave( 1,1) = -a1*zz
        wave(mu,1) = a1
        wave(mv,1) = 0.d0
        s(1) = -cc
    
        wave( 1,2) = a2*zz
        wave(mu,2) = a2
        wave(mv,2) = 0.d0
        s(2) = cc
        do m = 1,NEQNS
            amdq(m,i,iedge) = s(1)*wave(m,1)
            apdq(m,i,iedge) = s(2)*wave(m,2)
        enddo
    endif
    return
end subroutine rpn2_all_edges


! Note that valbig here is in SoA format
#ifdef CUDA
attributes(global) &
subroutine qad_gpu(valbig,mitot,mjtot,mbc,nvar, &
        svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy,&
        delt,mptr,max1d,cc,zz)

    implicit none

    integer, value, intent(in) :: mitot, mjtot, mbc, nvar, lenbc, max1d
    integer, value, intent(in) :: lratiox, lratioy, mptr
    double precision, value, intent(in) :: cc, zz
    double precision, value, intent(in) :: hx, hy, delt
    double precision, intent(in) :: valbig(mitot,mjtot,nvar)
    double precision, intent(in) :: qc1d(nvar,lenbc)
    double precision, intent(inout) :: svdflx(nvar,lenbc)

    double precision :: ql(nvar,max1d,SPACEDIM*2), qr(nvar,max1d,SPACEDIM*2)
    double precision :: amdq(nvar,max1d,SPACEDIM*2),  apdq(nvar,max1d,SPACEDIM*2)

    integer :: nc, nr
    integer :: ivar, istat
    integer :: i,j,ic,jc,l
    integer :: base
    integer :: tid
    double precision :: add

#ifdef DEBUG
    if (blockDim%y /= 1 .or. gridDim%y /= 1) then
        print *, "qad_gpu kernel should be called with blockDim%y == 1 and gridDim%y == 1"
        stop
    endif
#endif

    tid = (blockIdx%x-1) * blockDim%x + threadIdx%x
    nc = mjtot-2*mbc
    nr = mitot-2*mbc

    if (tid > (nc+nr)*2) then
        return
    endif
    ! prepare ql and qr
    !
    ! side 1
    base = 0
    if (tid > base .and. tid <= (base+nc)) then
        j = tid - base
        do ivar = 1, nvar
            ql(ivar,j,1) = valbig(mbc,j+mbc,ivar)
        enddo

        jc = (j-1)/lratioy+1
        do ivar = 1, nvar
            qr(ivar,j,1) = qc1d(ivar,jc)
        enddo
    endif

    ! side 2
    base = nc
    if (tid > base .and. tid <= (base+nr)) then
        i = tid - base
        do ivar = 1, nvar
            qr(ivar,i,2) = valbig(i+mbc,mjtot-mbc+1,ivar)
        enddo

        ic = (i-1)/lratiox+1 + nc/lratioy
        do ivar = 1, nvar
            ql(ivar,i,2) = qc1d(ivar,ic)
        enddo
    endif

    ! side 3
    base = nc + nr
    if (tid > base .and. tid <= (base+nc)) then
        j = tid - base
        do ivar = 1, nvar
            qr(ivar,j,3) = valbig(mitot-mbc+1,j+mbc,ivar)
        enddo

        jc = (j-1)/lratioy+1 + nc/lratioy + nr/lratiox
        do ivar = 1, nvar
            ql(ivar,j,3) = qc1d(ivar,jc)
        enddo
    endif

    ! side 4
    base = 2*nc + nr
    if (tid > base .and. tid <= (base+nr)) then
        i = tid - base
        do ivar = 1, nvar
            ql(ivar,i,4) = valbig(i+mbc,mbc,ivar)
        enddo

        ic = (i-1)/lratiox+1 + 2*nc/lratioy + nr/lratiox
        do ivar = 1, nvar
            qr(ivar,i,4) = qc1d(ivar,ic)
        enddo
    endif

    ! We only need amdq and apdq from this
    call rpn2_all_edges_dev(max1d,mbc, &
            ql,qr,amdq,apdq,tid,nc,nr,cc,zz)

    ! Need to make sure all tasks above are finished before
    ! Starting below
    ! But we can't synchronize arocss CUDA blocks (before CUDA 9.0)
    ! So we must rely on the execution order of CUDA threads
    !--------
    !  side 1
    !--------
    ! we have the wave. for side 1 add into sdflxm
    !
    base = 0
    if (tid > base .and. tid <= (base+nc)) then
        j = tid-base
        jc = (j-1)/lratioy + 1
        do ivar = 1, nvar
            add = amdq(ivar,j,1) * hy * delt &
                + apdq(ivar,j,1) * hy * delt
            istat = atomicadd(svdflx(ivar,jc), add)
        enddo
    endif

    !--------
    !  side 2
    !--------
    ! we have the wave. for side 2. add into sdflxp
    !
    base = nc
    if (tid > base .and. tid <= (base+nr)) then
        i = tid-base
        ic = (i-1)/lratiox + 1 + nc/lratioy
        do ivar = 1, nvar
            add = -amdq(ivar,i,2) * hx * delt &
                - apdq(ivar,i,2) * hx * delt
            istat = atomicadd(svdflx(ivar,ic), add)
        enddo
    endif

    !--------
    !  side 3
    !--------
    ! we have the wave. for side 3 add into sdflxp
    !
    base = nc + nr
    if (tid > base .and. tid <= (base+nc)) then
        j = tid-base
        jc = (j-1)/lratioy + 1 + nc/lratioy + nr/lratiox
        do ivar = 1, nvar
            add = - amdq(ivar,j,3) * hy * delt &
                - apdq(ivar,j,3) * hy * delt
            istat = atomicadd(svdflx(ivar,jc), add)
        enddo
    endif

    !--------
    !  side 4
    !--------
    ! we have the wave. for side 4. add into sdflxm
    !
    base = 2*nc + nr
    if (tid > base .and. tid <= (base+nr)) then
        i = tid-base
        ic = (i-1)/lratiox + 1 + 2*nc/lratioy + nr/lratiox
        do ivar = 1, nvar
            add = amdq(ivar,i,4) * hx * delt &
                + apdq(ivar,i,4) * hx * delt
            istat = atomicadd(svdflx(ivar,ic), add)
        enddo
    endif
    return
end subroutine qad_gpu

attributes(device) &
subroutine rpn2_all_edges_dev(maxm,mbc,ql,qr,amdq,apdq,tid,nc,nr,cc,zz)

    use amr_module

    implicit none

    integer, value, intent(in) :: maxm, mbc, tid, nc, nr
    double precision, value, intent(in) :: cc, zz
    double precision, intent(in)  ::   ql(NEQNS, 1:maxm,SPACEDIM*2)
    double precision, intent(in)  ::   qr(NEQNS, 1:maxm,SPACEDIM*2)
    double precision, intent(out) :: apdq(NEQNS, 1:maxm,SPACEDIM*2)
    double precision, intent(out) :: amdq(NEQNS, 1:maxm,SPACEDIM*2)

    double precision :: wave(NEQNS, NWAVES)
    double precision ::    s(1:NWAVES)

!     local arrays
!     ------------
    double precision :: delta1, delta2, delta3, a1, a2
    integer :: mu, mv, i, m, iedge, base

    ! side 1
    base = 0
    if (tid > base .and. tid <= (base+nc)) then
        iedge = 1
        mu = 2
        mv = 3
        i = tid - base
        delta1 = ql( 1,i,iedge) - qr( 1,i,iedge)
        delta2 = ql(mu,i,iedge) - qr(mu,i,iedge)
        a1 = (-delta1 + zz*delta2) / (2.d0*zz)
        a2 = (delta1 + zz*delta2) / (2.d0*zz)
    
    !        # Compute the waves.
    
        wave( 1,1) = -a1*zz
        wave(mu,1) = a1
        wave(mv,1) = 0.d0
        s(1) = -cc
    
        wave( 1,2) = a2*zz
        wave(mu,2) = a2
        wave(mv,2) = 0.d0
        s(2) = cc
        do m = 1,NEQNS
            amdq(m,i,iedge) = s(1)*wave(m,1)
            apdq(m,i,iedge) = s(2)*wave(m,2)
        enddo
    endif

    ! side 2
    base = nc
    if (tid > base .and. tid <= (base+nr)) then
        iedge = 2
        mu = 3
        mv = 2
        i = tid - base
        delta1 = ql( 1,i,iedge) - qr( 1,i,iedge)
        delta2 = ql(mu,i,iedge) - qr(mu,i,iedge)
        a1 = (-delta1 + zz*delta2) / (2.d0*zz)
        a2 = (delta1 + zz*delta2) / (2.d0*zz)
    
    !        # Compute the waves.
    
        wave( 1,1) = -a1*zz
        wave(mu,1) = a1
        wave(mv,1) = 0.d0
        s(1) = -cc
    
        wave( 1,2) = a2*zz
        wave(mu,2) = a2
        wave(mv,2) = 0.d0
        s(2) = cc
        do m = 1,NEQNS
            amdq(m,i,iedge) = s(1)*wave(m,1)
            apdq(m,i,iedge) = s(2)*wave(m,2)
        enddo
    endif

    ! side 3
    base = nc + nr
    if (tid > base .and. tid <= (base+nc)) then
        iedge = 3
        mu = 2
        mv = 3
        i = tid - base
        delta1 = ql( 1,i,iedge) - qr( 1,i,iedge)
        delta2 = ql(mu,i,iedge) - qr(mu,i,iedge)
        a1 = (-delta1 + zz*delta2) / (2.d0*zz)
        a2 = (delta1 + zz*delta2) / (2.d0*zz)
    
    !        # Compute the waves.
    
        wave( 1,1) = -a1*zz
        wave(mu,1) = a1
        wave(mv,1) = 0.d0
        s(1) = -cc
    
        wave( 1,2) = a2*zz
        wave(mu,2) = a2
        wave(mv,2) = 0.d0
        s(2) = cc
        do m = 1,NEQNS
            amdq(m,i,iedge) = s(1)*wave(m,1)
            apdq(m,i,iedge) = s(2)*wave(m,2)
        enddo
    endif

    ! side 4
    base = 2*nc + nr
    if (tid > base .and. tid <= (base+nr)) then
        iedge = 4
        mu = 3
        mv = 2
        i = tid - base
        delta1 = ql( 1,i,iedge) - qr( 1,i,iedge)
        delta2 = ql(mu,i,iedge) - qr(mu,i,iedge)
        a1 = (-delta1 + zz*delta2) / (2.d0*zz)
        a2 = (delta1 + zz*delta2) / (2.d0*zz)
    
    !        # Compute the waves.
    
        wave( 1,1) = -a1*zz
        wave(mu,1) = a1
        wave(mv,1) = 0.d0
        s(1) = -cc
    
        wave( 1,2) = a2*zz
        wave(mu,2) = a2
        wave(mv,2) = 0.d0
        s(2) = cc
        do m = 1,NEQNS
            amdq(m,i,iedge) = s(1)*wave(m,1)
            apdq(m,i,iedge) = s(2)*wave(m,2)
        enddo
    endif
    return
end subroutine rpn2_all_edges_dev
#endif


end module reflux_module
