! ==============================================================================
!  Regions Module
!   Module containing data structures and setup routines for region refinement.
!
! ==============================================================================
module regions_module

    implicit none
    save

    ! poly Region type definition
    type polyregion_type
        integer :: min_level,max_level
        integer :: num_points
        real(kind=8) :: t_low,t_hi
        real(kind=8) :: x_low,x_hi,y_low,y_hi
        real(kind=8), allocatable :: xpoints(:)
        real(kind=8), allocatable :: ypoints(:)

    end type polyregion_type


    integer :: num_regions
    type(polyregion_type), allocatable :: regions(:)
      
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
        write(parmunit,*) 'REGIONS:'
        write(parmunit,*) '-----------'

        if (present(fname)) then
            call opendatafile(unit,fname)
        else
            call opendatafile(unit,'regions.data')
        endif

        read(unit,"(i2)") num_regions
        if (num_regions == 0) then
            write(parmunit,*) '  No regions specified for refinement'
            
        else
            ! Refinement poly region data
            allocate(regions(num_regions))
            do i=1,num_regions
                read(unit,*) regions(i)%num_points, &
                             regions(i)%min_level, regions(i)%max_level, &
                             regions(i)%t_low, regions(i)%t_hi

                np = regions(i)%num_points
                allocate(regions(i)%xpoints(np+1))
                allocate(regions(i)%ypoints(np+1))

                do ip = 1,np
                    read(unit,*) regions(i)%xpoints(ip), &
                                 regions(i)%ypoints(ip) 
                enddo
                regions(i)%xpoints(np+1) = regions(i)%xpoints(1) 
                regions(i)%ypoints(np+1) = regions(i)%ypoints(1) 
                regions(i)%x_low = minval(regions(i)%xpoints)
                regions(i)%y_low = minval(regions(i)%ypoints)
                regions(i)%x_hi =  maxval(regions(i)%xpoints)
                regions(i)%y_hi =  maxval(regions(i)%ypoints)
            enddo
        endif
        close(unit)

    end subroutine set_regions


end module regions_module
