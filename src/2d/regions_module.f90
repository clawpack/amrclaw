! ==============================================================================
!  Regions Module
!   Module containing data structures and setup routines for region refinement.
!
! ==============================================================================
module regions_module

    implicit none
    save

    ! Region type definition
    type region_type
        integer :: min_level,max_level
        real(kind=8) :: x_low,y_low,x_hi,y_hi,t_low,t_hi
    end type region_type

    ! poly Region type definition
    type polyregion_type
        integer :: min_level,max_level
        integer :: num_points
        real(kind=8) :: t_low,t_hi
        real(kind=8), allocatable :: xpoints(:)
        real(kind=8), allocatable :: ypoints(:)

    end type polyregion_type

    integer :: num_regions
    type(region_type), allocatable :: regions(:)


    integer :: num_polyregions
    type(polyregion_type), allocatable :: polyregions(:)
      
contains


    subroutine set_regions(fname)
        ! this sets arbitrary polygon regions

        use amr_module
      
        implicit none
      
        ! Function Arguments
        character(len=*), optional, intent(in) :: fname
      
        ! Locals
        integer, parameter :: unit = 7
        integer :: i,ip,np

        write(parmunit,*) ' '
        write(parmunit,*) '--------------------------------------------'
        write(parmunit,*) 'POLY REGIONS:'
        write(parmunit,*) '-----------'

        if (present(fname)) then
            call opendatafile(unit,fname)
        else
            call opendatafile(unit,'regions.data')
        endif

        read(unit,"(i2)") num_polyregions
        if (num_polyregions == 0) then
            write(parmunit,*) '  No polyregions specified for refinement'
            
        else
            ! Refinement poly region data
            allocate(polyregions(num_polyregions))
            do i=1,num_polyregions
                read(unit,*) polyregions(i)%num_points, &
                             polyregions(i)%min_level, polyregions(i)%max_level, &
                             polyregions(i)%t_low, polyregions(i)%t_hi

                np = polyregions(i)%num_points
                allocate(polyregions(i)%xpoints(np))
                allocate(polyregions(i)%ypoints(np))
                do ip = 1,np
                    read(unit,*) polyregions(i)%xpoints(ip), &
                                 polyregions(i)%ypoints(ip) 
                enddo
            enddo
        endif
        close(unit)

    end subroutine set_regions


end module regions_module
