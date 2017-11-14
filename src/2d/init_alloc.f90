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
    integer :: istat
#endif
    
!    if (.not.allocated(storage)) then   ! old way, changed mjb sept. 2014
    if (.not.allocated(alloc)) then      ! new way, use allocatable arrays, not pointers
        memsize = 1000000
!        allocate(storage(memsize)) 
#ifdef CUDA
        allocate(alloc(memsize), stat=istat, pinned=plog)
        if (.not. plog) then
            print *, "Warning: allocating pinned memory in init_alloc() failed"
        endif
#else
        allocate(alloc(memsize))
#endif
!        alloc => storage
        print *, "Storage allocated..."
    else
        print *, "Storage already allocated!"
    endif

end subroutine init_alloc

