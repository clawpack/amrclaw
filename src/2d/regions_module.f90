! ==============================================================================
!  Regions Module
!   Module containing data structures and setup routines for region refinement.
!
! ==============================================================================
module regions_module

    implicit none
    save

    ! Region type definition (old style rectangles still supported for now)
    type region_type
        integer :: min_level,max_level
        real(kind=8) :: x_low,y_low,x_hi,y_hi,t_low,t_hi
    end type region_type

    ! New flag regions specified as Ruled Rectangles
    type ruled_region_type
        integer :: min_level,max_level, ixy, method, nrules
        real(kind=8) :: t_low,t_hi,ds
        real(kind=8) :: x1bb,x2bb,y1bb,y2bb  ! bounding box
        real(kind=8), allocatable, dimension(:) :: s, lower, upper
        character(len=200) :: file_name, name
    end type ruled_region_type

    logical, private :: module_setup

    ! old style regions:
    integer :: num_regions
    type(region_type), allocatable  :: regions(:)
      
    ! new style flagregions
    integer :: num_rregions
    type(ruled_region_type), target, allocatable :: rregions(:)
      
contains

    subroutine set_regions(fname)

        ! Read old-style rectangular regions from regions.data
        
        use amr_module
      
        implicit none
      
        ! Function Arguments
        character(len=*), optional, intent(in) :: fname
      
        ! Locals
        integer, parameter :: unit = 7
        integer :: i

        if (.not. module_setup) then

            write(parmunit,*) ' '
            write(parmunit,*) '--------------------------------------------'
            write(parmunit,*) 'REGIONS:'
            write(parmunit,*) '-----------'

            if (present(fname)) then
                call opendatafile(unit,fname)
            else
                call opendatafile(unit,'regions.data')
            endif

            read(unit,*) num_regions
            if (num_regions == 0) then
                write(parmunit,*) '  No regions specified for refinement'
                
            else
                ! Refinement region data
                allocate(regions(num_regions))
                do i=1,num_regions
                    read(unit,*) regions(i)%min_level, regions(i)%max_level, &
                                 regions(i)%t_low, regions(i)%t_hi, &
                                 regions(i)%x_low, regions(i)%x_hi, &
                                 regions(i)%y_low, regions(i)%y_hi
                enddo
            endif
            close(unit)

            call set_rregions()

            module_setup = .true.
        end if

    end subroutine set_regions


    subroutine set_rregions(fname)

        ! Read new-style Ruled Rectangles from flagregions.data
        
        use amr_module
      
        implicit none
      
        ! Function Arguments
        character(len=*), optional, intent(in) :: fname
      
        ! Locals
        integer, parameter :: unit = 7
        integer, parameter :: unit2 = 45
        integer :: i,j,nrules1, spatial_region_type
        logical :: foundFile
        type(ruled_region_type), pointer :: rr
        real(kind=8) :: rr_x1,rr_x2,rr_y1,rr_y2


        write(parmunit,*) ' '
        write(parmunit,*) '--------------------------------------------'
        write(parmunit,*) 'RULED REGIONS:'
        write(parmunit,*) '-----------'

        if (present(fname)) then
            call opendatafile(unit,fname)
        else
            call opendatafile(unit,'flagregions.data')
        endif

        read(unit,*) num_rregions
        if (num_rregions == 0) then
            write(parmunit,*) '  No ruled regions specified for refinement'
            
        else
            ! Refinement region data
            allocate(rregions(num_rregions))
            do i=1,num_rregions
                rr => rregions(i)
                read(unit,*) rr%name
                rr%name = trim(rr%name)
                read(unit,*) rr%min_level
                read(unit,*) rr%max_level
                read(unit,*) rr%t_low
                read(unit,*) rr%t_hi
                read(unit,*) spatial_region_type

                if (spatial_region_type == 1) then
                    ! read rectangle extent:
                    read(unit,*) rr_x1,rr_x2,rr_y1,rr_y2
                    ! turn into ruled rectangle:
                    rr%ixy = 1
                    rr%method = 0
                    rr%nrules = 2
                    allocate(rr%s(rr%nrules), rr%lower(rr%nrules), & 
                             rr%upper(rr%nrules))
                    rr%s(1) = rr_x1
                    rr%s(2) = rr_x2
                    rr%ds = rr_x2 - rr_x1
                    rr%lower(1) = rr_y1
                    rr%upper(1) = rr_y2
                    rr%lower(2) = rr_y1
                    rr%upper(2) = rr_y2
                
                else if (spatial_region_type == 2) then    
                    read(unit,*) rr%file_name
                    write(6,*) '+++ Ruled region name: ', rr%name
                    write(6,*) '+++ Ruled region file_name: ', rr%file_name
                    inquire(file=trim(rr%file_name),exist=foundFile)
                    if (.not. foundFile) then
                      write(*,*) 'Missing rregions file...'
                      write(*,*) 'Looking for: ',trim(rr%file_name)
                      stop
                      endif

                    open(unit=unit2,file=trim(rr%file_name),status='old')
                    read(unit2,*) rr%ixy
                    read(unit2,*) rr%method
                    read(unit2,*) rr%ds
                    read(unit2,*) rr%nrules
                    allocate(rr%s(rr%nrules), rr%lower(rr%nrules), & 
                             rr%upper(rr%nrules))
                    do j=1,rr%nrules
                        read(unit2,*) rr%s(j), rr%lower(j), rr%upper(j)
                        enddo
                        
                else
                    write(6,*) '*** Error: unexpected spatial_region_type'
                endif
                    
                ! compute bounding box:
                
                if (rr%method == 0) then
                    nrules1 = rr%nrules - 1
                else
                    nrules1 = rr%nrules
                    endif
                    
                if (rr%ixy == 1) then
                    rr%x1bb = rr%s(1)
                    rr%x2bb = rr%s(rr%nrules)
                    rr%y1bb = minval(rr%lower(1:nrules1))
                    rr%y2bb = maxval(rr%upper(1:nrules1))
                else
                    rr%y1bb = rr%s(1)
                    rr%y2bb = rr%s(rr%nrules)
                    rr%x1bb = minval(rr%lower(1:nrules1))
                    rr%x2bb = maxval(rr%upper(1:nrules1))
                    endif
                    
                write(6,*) '+++ rregion bounding box: '
                write(6,*) rr%x1bb,rr%x2bb,rr%y1bb,rr%y2bb
                
                write(6,*) '+++ i, rr%s(1), rr%ds: ',i, rr%s(1), rr%ds
                enddo
                
            endif
            close(unit)


    end subroutine set_rregions

end module regions_module
