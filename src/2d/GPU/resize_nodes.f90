! ============================================================================
!  Program:     AMRClaw
!  File:        resize_nodes.f90
!  Author:      Mandli and Marsha Berger
! ============================================================================
!  Description:  Resize the node space if have run out 
! ============================================================================

subroutine resize_nodes(new_size,status)
    
    use amr_module
    implicit none
    
    integer, intent(out) :: status
    integer, intent(in) :: new_size

    integer :: i
    
    real(CLAW_REAL), allocatable, target, dimension(:,:) :: new_rnode
    integer, allocatable, target, dimension(:,:) :: new_node
    integer, allocatable, target, dimension(:) :: new_listOfGrids

    type(gpu_3d_real_ptr_type), allocatable, target, dimension(:) :: new_grid_data_d
    type(gpu_3d_real_ptr_type), allocatable, target, dimension(:) :: new_grid_data_d_copy2
    type(gpu_3d_real_ptr_type), allocatable, target, dimension(:) :: new_aux_d

    type(cpu_2d_int_ptr_type), allocatable, target, dimension(:) :: new_cflux_hh
    type(gpu_2d_int_ptr_type), allocatable, target, dimension(:) :: new_cflux_hd
    type(cpu_1d_real_ptr_type), allocatable, target, dimension(:) :: new_fflux_hh
    type(gpu_1d_real_ptr_type), allocatable, target, dimension(:) :: new_fflux_hd

    print *, "Expanding maximum number of grids from ", maxgr," to ", new_size

    ! first for rnode
    allocate(new_rnode(rsize,new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_rnode(1:rsize,1:maxgr) = rnode     ! new way, use allocatable, not pointer       

    call move_alloc(new_rnode,rnode)

    ! next for node
    allocate(new_node(nsize,new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_node(1:nsize,1:maxgr) = node     ! new way, use allocatable, not pointer       

    call move_alloc(new_node,node)

    !! need to rethread new space to be able to use it when new grids requested
    do i = maxgr+1, new_size
       node(nextfree,i) = i+1
    end do
    ! reset last one to null
    node(nextfree, new_size) = clawpack_null

    ! next for listOfGrids
    allocate(new_listOfGrids(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_listOfGrids(1:maxgr) = listOfGrids     ! new way, use allocatable, not pointer       

    call move_alloc(new_listOfGrids,listOfGrids)

#ifdef CUDA
    allocate(new_grid_data_d(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_grid_data_d(1:maxgr) = grid_data_d
    call move_alloc(new_grid_data_d, grid_data_d)


    allocate(new_grid_data_d_copy2(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_grid_data_d_copy2(1:maxgr) = grid_data_d_copy2
    call move_alloc(new_grid_data_d_copy2, grid_data_d_copy2)

    allocate(new_aux_d(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_aux_d(1:maxgr) = aux_d
    call move_alloc(new_aux_d, aux_d)

    deallocate(waveSpeedsX)
    allocate(waveSpeedsX(ws_len,new_size),STAT=status)
    if (status > 0) then
        return
    endif

    deallocate(waveSpeedsY)
    allocate(waveSpeedsY(ws_len,new_size),STAT=status)
    if (status > 0) then
        return
    endif

    allocate(new_cflux_hh(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_cflux_hh(1:maxgr) = cflux_hh
    call move_alloc(new_cflux_hh, cflux_hh)

    allocate(new_cflux_hd(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_cflux_hd(1:maxgr) = cflux_hd
    call move_alloc(new_cflux_hd, cflux_hd)

    deallocate(cflux_dd)
    allocate(cflux_dd(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    cflux_dd = cflux_hd

    allocate(new_fflux_hh(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_fflux_hh(1:maxgr) = fflux_hh
    call move_alloc(new_fflux_hh, fflux_hh)

    allocate(new_fflux_hd(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    new_fflux_hd(1:maxgr) = fflux_hd
    call move_alloc(new_fflux_hd, fflux_hd)

    deallocate(fflux_dd)
    allocate(fflux_dd(new_size),STAT=status)
    if (status > 0) then
        return
    endif
    fflux_dd = fflux_hd
#endif

    ! reset maxgr and next free node,  to continue
    ndfree = maxgr + 1     
    
    maxgr = new_size
    return
    
end subroutine resize_nodes
