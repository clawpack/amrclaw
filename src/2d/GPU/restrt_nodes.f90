! ============================================================================
!  Program:     AMRClaw
!  File:        restrt_alloc.f90
!  Created:     2009-10-22
!  Author:      Marsha Berger and Randy LeVeque
! ============================================================================
!  Description:  Initialization of alloc storage for restart
! ============================================================================


subroutine restrt_nodes(isize)

    use amr_module
    implicit none
    integer :: isize

    maxgr = isize

    if (.not.allocated(rnode)) then     ! allocatable nodal arrays now 
        write(6,*)"allocating ",maxgr," -sized rnode array"
        allocate(rnode(rsize,maxgr))
        print *, "rnode storage allocated of size ",maxgr," at restart"
    else
        print *, "rnode storage already allocated!"
    endif

    if (.not.allocated(node)) then     ! allocatable nodal arrays now 
        write(6,*)"allocating ",maxgr," -sized node array"
        allocate(node(nsize,maxgr))
        print *, "node storage allocated of size ",maxgr," at restart"
    else
        print *, "node storage already allocated!"
    endif

    if (.not.allocated(listOfGrids)) then     ! allocatable nodal arrays now 
        write(6,*)"allocating ",maxgr," -sized listOfGrids "
        allocate(listOfGrids(maxgr))
        print *, "listOfGrids allocated of size ",maxgr," at restart"
    else
        print *, "listOfGrids already allocated!"
    endif

#ifdef CUDA
    if (.not.allocated(grid_data_d)) then   
        write(6,*)"allocating ",maxgr," -sized grid_data_d"
        allocate(grid_data_d(maxgr))
        print *, "grid_data_d storage allocated of size ",maxgr," at restart"
    else
        print *, "grid_data_d already allocated!"
    endif

    if (.not.allocated(grid_data_d_copy2)) then   
        write(6,*)"allocating ",maxgr," -sized grid_data_d_copy2"
        allocate(grid_data_d_copy2(maxgr))
        print *, "grid_data_d_copy2 storage allocated of size ",maxgr," at restart"
    else
        print *, "grid_data_d_copy2 already allocated!"
    endif

    if (.not.allocated(aux_d)) then   
        write(6,*)"allocating ",maxgr," -sized aux_d"
        allocate(aux_d(maxgr))
        print *, "aux_d storage allocated of size ",maxgr," at restart"
    else
        print *, "aux_d already allocated!"
    endif

    if (.not.allocated(waveSpeedsX)) then   
        write(6,*)"allocating ",maxgr," -sized waveSpeedsX"
        allocate(waveSpeedsX(ws_len,maxgr))
        print *, "waveSpeedsX storage allocated of size ",maxgr," at restart"
    else
        print *, "waveSpeedsX already allocated!"
    endif

    if (.not.allocated(waveSpeedsY)) then   
        write(6,*)"allocating ",maxgr," -sized waveSpeedsY"
        allocate(waveSpeedsY(ws_len,maxgr))
        print *, "waveSpeedsY storage allocated of size ",maxgr," at restart"
    else
        print *, "waveSpeedsY already allocated!"
    endif

    print *, "allocating fflux and cflux at restart ..."
    allocate(fflux_hh(maxgr))
    allocate(fflux_hd(maxgr))
    allocate(fflux_dd(maxgr))
    allocate(cflux_hh(maxgr))
    allocate(cflux_hd(maxgr))
    allocate(cflux_dd(maxgr))
#endif

end subroutine restrt_nodes


