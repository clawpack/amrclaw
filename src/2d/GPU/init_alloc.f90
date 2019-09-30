! ============================================================================
!  Program:     AMRClaw
!  File:        init_alloc.f90
!  Created:     2009-01-21
!  Author:      Kyle Mandli and Marsha Berger
! ============================================================================
!  Description:  Initialization of alloc storage
! ============================================================================


subroutine init_alloc()
    
    use amr_module
    implicit none

#ifdef CUDA
    logical :: plog
#endif
    
    if (.not.allocated(alloc)) then      ! new way, use allocatable arrays, not pointers
        memsize = 4000000
#ifdef CUDA
        allocate(alloc(memsize), pinned=plog)
        if (.not. plog) then
            print *, "Warning: allocating pinned memory in init_alloc() failed"
        endif
#else
!        allocate(storage(memsize)) 
        allocate(alloc(memsize))
!        alloc => storage
#endif
        print *, "Storage allocated..."
    else
        print *, "Storage already allocated!"
    endif
    
end subroutine init_alloc

