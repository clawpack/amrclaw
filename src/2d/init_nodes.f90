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
    
    
end subroutine init_nodes

