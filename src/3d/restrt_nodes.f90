! ============================================================================
!  Program:     AMRClaw
!  File:        restrt_nodes.f90
!  Created:     2018-06-15
!  Author:      Marsha Berger
! ============================================================================
!  Description:  Initialization of nodal storage for restart
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

end subroutine restrt_nodes


