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
    
    maxgr = 50000   
    if (.not.allocated(rnode)) then      ! new way, use allocatable arrays, not pointers
        write(*,*) "rsize ",rsize
        allocate(rnode(rsize,maxgr))
        print *, "rnode allocated..."
    else
        print *, "rnode already allocated!"
    endif
    if (.not.allocated(node)) then      ! new way, use allocatable arrays, not pointers
        write(*,*) "nsize ", nsize
        allocate(node(nsize,maxgr))
        print *, "node allocated..."
    else
        print *, "rnode already allocated!"
    endif
    
    
end subroutine init_nodes

