! ============================================================================
!  Program:     AMRClaw
!  File:        restrt_alloc.f90
!  Created:     2009-10-22
!  Author:      Marsha Berger and Randy LeVeque
! ============================================================================
!  Description:  Initialization of alloc storage for restart
! ============================================================================


subroutine restrt_alloc(isize)

    use amr_module
    implicit none
    integer :: isize


!    if (.not.allocated(storage)) then  ! old way, changed mjb Sept. 2014
    if (.not.allocated(alloc)) then     ! new way uses allocatable, not pointer, arrays
        write(6,*)"allocating ",isize," -sized alloc array"
        memsize = isize
!       allocate(storage(memsize))
        allocate(alloc(memsize))
!       alloc => storage
        print *, "Storage allocated of size ",memsize," at restart"
    else
        print *, "Storage already allocated!"
    endif

end subroutine restrt_alloc


