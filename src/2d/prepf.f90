!
!> For new fine grids (grids on level **level**), allocate space for 
!! saving fluxes at boundary (even if it's physical boundary) after 
!! each integration step.
!! These fluxes are for future conservative fixing of coarse grids
!! (level **level**-1 grids).
!! The address of this space is stored in node(ffluxptr, mkid) for grid
!! mkid. 
!!
!! Note that if the refinment ratio is r, fluxes from every r cells 
!! on grid **mkid** are saved to one slot in this space, since every r
!! cells on grid **mkid** border one cell on level **level**-1 grids.
!! \param[in] level boudnary lists of grids on this level get updated
!! \param[in] nvar number of equations for the system
!! \param[in] naux number of auxiliary variables for the system
! ----------------------------------------------------------
!
#include "amr_macros.H"
subroutine prepf(level,nvar,naux)
    !
    use amr_module
#ifdef CUDA
    use memory_module, only: gpu_allocate
    use cuda_module, only: device_id
#endif
    implicit double precision (a-h,o-z)

    !
    ! ::::::::::::::::::::::::::: PREPF :::::::::::::::::::::::::::::
    !
    ! prepare new fine grids to save fluxes after each integration step
    ! for future conservative fixing of coarse grids
    ! save all boundary fluxes of fine grid (even if phys. bndry.) -
    ! but only save space for every intrat. (remember - 4 fluxes).
    !
    ! :::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::::::
    !
    mkid = lstart(level)
    do while (mkid .ne. 0)
        nx    = node(ndihi,mkid)-node(ndilo,mkid) + 1
        ny    = node(ndjhi,mkid)-node(ndjlo,mkid) + 1
        ikeep = nx/intratx(level-1)
        jkeep = ny/intraty(level-1)
        lenbc = 2*(ikeep+jkeep)
        !
        !         get twice the storage, one for plus or minus fluxes, the
        !         other for coarse solution for wave fixing. also for aux vars.
        !
        node(ffluxptr,mkid) = igetsp(2*nvar*lenbc+naux*lenbc)
        ist = node(ffluxptr,mkid)
#ifdef CUDA
        ! TODO: I should check whether zeroize the array is necessary here
        call gpu_allocate(node_data(mkid,FFLUXPTR_D)%dataptr, &
        device_id, 1, 2*nvar*lenbc+naux*lenbc)
        print *, "Allocate ffluxptr_d of grid: ", mkid
        print *, "at: ", loc(node_data(mkid,FFLUXPTR_D)%dataptr)
#endif

        do i = 1, 2*nvar*lenbc + naux*lenbc
            alloc(ist+i-1) = 0.d0
        enddo
        mkid = node(levelptr,mkid)
    enddo

    return
end subroutine prepf
