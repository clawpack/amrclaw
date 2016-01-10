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

!      if (.not.allocated(storage)) then   ! old way, changed mjb jan. 2015
       if (.not.allocated(alloc)) then     ! new way, use allocatable arrays, not pointers
          memsize = 4000000
!         allocate(storage(memsize))
          allocate(alloc(memsize))
!          alloc => storage
          print *, "Storage allocated..."
      else
          print *, "Storage already allocated!"
      endif
    
      end subroutine init_alloc

