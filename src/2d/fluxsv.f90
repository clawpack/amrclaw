!> When a coarse grid cell is advanced, if it borders a coarse-fine
!! interface, the flux or wave that emits from the interface and goes into the coarse cell 
!! is added to the corresponding location in **node(ffluxptr, mkid)**
!! for conservative fix later, where **mkid** is grid number of the
!! nested fine grid on the other side of that interface. 
!! 
!! 
! ----------------------------------------------------------
!
#include "amr_macros.H"

#ifndef CUDA
subroutine fluxsv(mptr,xfluxm,xfluxp,yfluxm,yfluxp,listbc, &
        ndimx,ndimy,nvar,maxsp,dtc,hx,hy)

    use amr_module
    implicit double precision (a-h,o-z)


    dimension xfluxp(nvar,ndimx,ndimy), yfluxp(nvar,ndimx,ndimy)
    dimension xfluxm(nvar,ndimx,ndimy), yfluxm(nvar,ndimx,ndimy)
    dimension listbc(5,maxsp)
    !
    ! :::::::::::::::::::: FLUXSV :::::::::::::::::::::::::
    !
    !  coarse grids should save their fluxes in cells adjacent to
    !  their nested fine grids, for later conservation fixing.
    !  listbc holds info for where to save which fluxes.
    !  xflux holds 'f' fluxes, yflux holds 'g' fluxes.
    !
    ! :::::::::::::::::::::::::::::;:::::::::::::::::::::::


    level   = node(nestlevel,mptr)

    do ispot = 1,maxsp
        if (listbc(1,ispot).eq.0) then 
            return
        else
            !
            mkid     = listbc(4,ispot)
            intopl   = listbc(5,ispot)
            nx       = node(ndihi,mkid) - node(ndilo,mkid) + 1
            ny       = node(ndjhi,mkid) - node(ndjlo,mkid) + 1
            kidlst   = node(ffluxptr,mkid)
            i        = listbc(1,ispot)
            j        = listbc(2,ispot)
            inlist   = kidlst + nvar*(intopl-1) - 1
            !
            ! side k (listbc 3) has which side of coarse cell has interface
            ! so can save appropriate fluxes.  (dont know why we didnt have
            ! which flux to save directly (i.e. put i+1,j to save that flux
            ! rather than putting in cell center coords).

            if (listbc(3,ispot) .eq. 1) then
                !           ::::: Cell i,j is on right side of a fine grid
                do ivar = 1, nvar
                    alloc(inlist + ivar) = -xfluxp(ivar,i,j)*dtc*hy
                enddo
                !         write(dbugunit,901) i,j,1,(xfluxp(ivar,i,j),ivar=1,nvar)
            endif

            if (listbc(3,ispot) .eq. 2) then
                !           ::::: Cell i,j on bottom side of fine grid
                do ivar = 1, nvar
                    alloc(inlist + ivar) = -yfluxm(ivar,i,j+1)*dtc*hx
                enddo
                !         write(dbugunit,901) i,j,2,(yfluxm(ivar,i,j+1),ivar=1,nvar)
            endif

            if (listbc(3,ispot) .eq. 3) then
                !           ::::: Cell i,j on left side of fine grid
                do ivar = 1, nvar
                    alloc(inlist + ivar) = -xfluxm(ivar,i+1,j)*dtc*hy
                enddo
                !         write(dbugunit,901) i,j,3,(xfluxm(ivar,i+1,j),ivar=1,nvar)
            endif

            if (listbc(3,ispot) .eq. 4) then
                !           ::::: Cell i,j on top side of fine grid
                do ivar = 1, nvar
                    alloc(inlist + ivar) = -yfluxp(ivar,i,j)*dtc*hx
                enddo
                !         write(dbugunit,901) i,j,4,(yfluxp(ivar,i,j),ivar=1,nvar)
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
                    alloc(inlist + ivar) = yfluxm(ivar,i,j+1)*dtc*hx
                enddo
                !         write(dbugunit,901) i,j,5,(yfluxm(ivar,i,j+1),ivar=1,nvar)
                !  901        format(2i4," side",i3,4e15.7)
            endif
            !
            if (listbc(3,ispot) .eq. 6) then
                !           ::::: Cell i,j on bottom side of fine grid with spherical mapped bc
                do ivar = 1, nvar
                    alloc(inlist + ivar) = yfluxp(ivar,i,j)*dtc*hx
                enddo
                !         write(dbugunit,901) i,j,6,(yfluxp(ivar,i,j),ivar=1,nvar)
            endif
        endif
    enddo
    return
end subroutine fluxsv
#endif

#ifdef CUDA

! This is the CPU version of fluxsv_gpu
! It's mainly for development of GPU version
subroutine fluxsv_cpu(mptr,&
                xfluxm,xfluxp,yfluxm,yfluxp,&
                listbc,&
                fflux,&
                ndimx,ndimy,nvar,maxsp,dtc,hx,hy)

    use amr_module
    implicit none


    integer, intent(in) :: mptr
    integer, intent(in) :: ndimx, ndimy, nvar, maxsp
    double precision, intent(in) :: dtc, hx, hy
    double precision, intent(in) :: xfluxp(ndimx,ndimy,nvar), yfluxp(ndimx,ndimy,nvar)
    double precision, intent(in) :: xfluxm(ndimx,ndimy,nvar), yfluxm(ndimx,ndimy,nvar)
    integer, intent(in) :: listbc(5,maxsp)
    type(cpu_1d_real_ptr_type), intent(in) :: fflux(15000)

    ! local
    integer :: ispot, mkid, intopl, loc
    integer :: i,j, ivar
    integer :: tid
    integer :: kidlst, inlist

    ! :::::::::::::::::::: FLUXSV :::::::::::::::::::::::::
    !
    !  coarse grids should save their fluxes in cells adjacent to
    !  their nested fine grids, for later conservation fixing.
    !  listbc holds info for where to save which fluxes.
    !  xflux holds 'f' fluxes, yflux holds 'g' fluxes.
    !
    ! :::::::::::::::::::::::::::::;:::::::::::::::::::::::

    do tid = 1,maxsp
        if (listbc(1,tid) .eq. 0) then 
            cycle
        endif

        ispot = tid
        
        mkid     = listbc(4,ispot)
        intopl   = listbc(5,ispot)
        i        = listbc(1,ispot)
        j        = listbc(2,ispot)

        loc = nvar*(intopl-1)

        ! side k (listbc 3) has which side of coarse cell has interface
        ! so can save appropriate fluxes.  (dont know why we didnt have
        ! which flux to save directly (i.e. put i+1,j to save that flux
        ! rather than putting in cell center coords).

        ! data in fflux(mkid)%flux_fix is in AoS format. Namely,
        ! f(ivar, i, j)
        if (listbc(3,ispot) .eq. 1) then
            !           ::::: Cell i,j is on right side of a fine grid
            do ivar = 1, nvar
                fflux(mkid)%ptr(loc + ivar) = -xfluxp(i,j,ivar)*dtc*hy
                ! alloc(inlist + ivar) = -xfluxp(i,j,ivar)*dtc*hy
            enddo
        endif

        if (listbc(3,ispot) .eq. 2) then
            !           ::::: Cell i,j on bottom side of fine grid
            do ivar = 1, nvar
                fflux(mkid)%ptr(loc + ivar) = -yfluxm(i,j+1,ivar)*dtc*hx
                ! alloc(inlist + ivar) = -yfluxm(i,j+1,ivar)*dtc*hx
            enddo
        endif

        if (listbc(3,ispot) .eq. 3) then
            !           ::::: Cell i,j on left side of fine grid
            do ivar = 1, nvar
                fflux(mkid)%ptr(loc + ivar) = -xfluxm(i+1,j,ivar)*dtc*hy
                ! alloc(inlist + ivar) = -xfluxm(i+1,j,ivar)*dtc*hy
            enddo
        endif

        if (listbc(3,ispot) .eq. 4) then
            !           ::::: Cell i,j on top side of fine grid
            do ivar = 1, nvar
                fflux(mkid)%ptr(loc + ivar) = -yfluxp(i,j,ivar)*dtc*hx 
                ! alloc(inlist + ivar) = -yfluxp(i,j,ivar)*dtc*hx
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
                ! alloc(inlist + ivar) = yfluxm(i,j+1,ivar)*dtc*hx
            enddo
        endif
        !
        if (listbc(3,ispot) .eq. 6) then
            !           ::::: Cell i,j on bottom side of fine grid with spherical mapped bc
            do ivar = 1, nvar
                fflux(mkid)%ptr(loc + ivar) = yfluxp(i,j,ivar)*dtc*hx
                ! alloc(inlist + ivar) = yfluxp(i,j,ivar)*dtc*hx
            enddo
        endif
    enddo
    return
end subroutine fluxsv_cpu

#endif
