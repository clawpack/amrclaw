! ============================================================================
!  Program:     AMRClaw
!  File:        init_nodes.f90
!  Created:     2009-01-21
!  Author:      Marsha Berger
! ============================================================================
!  Description:  Initialization of rnode and node storage
! ============================================================================


subroutine init_nodes()
    
    use amr_module
    implicit none
    
    maxgr = 10000   
    if (.not.allocated(rnode)) then      ! new way, use allocatable arrays, not pointers
        allocate(rnode(rsize,maxgr))
        print *, "rnode allocated..."
    else
        print *, "rnode already allocated!"
    endif
    if (.not.allocated(node)) then      ! new way, use allocatable arrays, not pointers
        allocate(node(nsize,maxgr))
        print *, "node allocated..."
    else
        print *, "rnode already allocated!"
    endif
    if (.not.allocated(listOfGrids)) then   ! include all vars whose size depends on maxgr
        allocate(listOfGrids(maxgr))
        print *, "listOfGrids allocated..."
    else
        print *, "listOfGrids already allocated!"
    endif
    
#ifdef CUDA
    if (.not.allocated(grid_data_d)) then   
        allocate(grid_data_d(maxgr))
        print *, "grid_data_d allocated..."
    else
        print *, "grid_data_d already allocated!"
    endif

    if (.not.allocated(grid_data_d_copy2)) then   
        allocate(grid_data_d_copy2(maxgr))
        print *, "grid_data_d_copy2 allocated..."
    else
        print *, "grid_data_d_copy2 already allocated!"
    endif

    if (.not.allocated(aux_d)) then   
        allocate(aux_d(maxgr))
        print *, "aux_d allocated..."
    else
        print *, "aux_d already allocated!"
    endif

    if (.not.allocated(waveSpeedsX)) then   
        allocate(waveSpeedsX(ws_len,maxgr))
        print *, "waveSpeedsX on device allocated..."
    else
        print *, "waveSpeedsX already allocated!"
    endif

    if (.not.allocated(waveSpeedsY)) then   
        allocate(waveSpeedsY(ws_len,maxgr))
        print *, "waveSpeedsY on device allocated..."
    else
        print *, "waveSpeedsY already allocated!"
    endif

    print *, "allocating fflux and cflux ..."
    allocate(fflux_hh(maxgr))
    allocate(fflux_hd(maxgr))
    allocate(fflux_dd(maxgr))
    allocate(cflux_hh(maxgr))
    allocate(cflux_hd(maxgr))
    allocate(cflux_dd(maxgr))
#endif
    
end subroutine init_nodes

