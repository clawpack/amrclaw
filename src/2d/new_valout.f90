!> Output the results for a general system of conservation laws
!! in 2 dimensions
!!
!! Write the results to the file fort.q<iframe>
!! Use format required by matlab script  plotclaw2.m or Python tools
!!
!! set outaux = .true. to also output the aux arrays to fort.a<iframe>
subroutine new_valout(level_begin, level_end, time, num_eqn, num_aux)

    use amr_module, only: t0, output_aux_onlyonce, output_aux_components
    use amr_module, only: frame => matlabu, num_grids => ngrids
    use amr_module, only: num_ghost => nghost, lstart, output_format
    use amr_module, only: hxposs, hyposs,

#ifdef NETCDF
    use netcdf
#endif

    implicit none

    ! Input
    integer :: level_begin, level_end, num_eqn, num_aux
    real(kind=8) :: time

    ! Locals
    integer, parameter :: out_unit = 50
    integer :: i, j, m, level, output_aux_num, num_stop, digit, num_dim
    integer :: level_ptr, num_cells(2), num_grids
    real(kind=8) :: lower_corner(2), delta(2)
    logical :: out_aux
    character(len=10) :: file_name(5)

    integer :: clock_start, clock_finish, clock_rate
    real(kind=8) cpu_start, cpu_finish

    character(len=*), parameter :: q_header_format_1d =                        &
                                    "(i6,'                 grid_number',/," // &
                                     "i6,'                 AMR_level',/,"   // &
                                     "i6,'                 mx',/,"          // &
                                     "e26.16,'    xlow', /, "               // &
                                     "e26.16,'    dx',/)"
    character(len=*), parameter :: q_header_format_2d =                        &
                                    "(i6,'                 grid_number',/," // &
                                     "i6,'                 AMR_level',/,"   // &
                                     "i6,'                 mx',/,"          // &
                                     "i6,'                 my',/"           // &
                                     "e26.16,'    xlow', /, "               // &
                                     "e26.16,'    ylow', /,"                // &
                                     "e26.16,'    dx', /,"                  // &
                                     "e26.16,'    dy',/)"
    character(len=*), parameter :: q_file_format = "()"
    character(len=*), parameter :: t_file_format = "(e18.8,'    time', /,"  // &
                                           "i6,'                 meqn'/,"   // &
                                           "i6,'                 ngrids'/," // &
                                           "i6,'                 naux'/,"   // &
                                           "i6,'                 ndim'/,"   // &
                                           "i6,'                 nghost'/,/)"
    character(len=*), parameter :: console_format = "('AMRCLAW: Frame " // &
                              "',i4,' output files done at time t = ', d13.6,/)"

    !
    !
    pure integer function iadd(ivar, i, j)
        implicit none
        iadd = q_loc + ivar - 1 + nvar*((j-1)*mitot+i-1)
    end function iadd

    !
    !
    pure integer function iaddaux(ivar, i, j)
        implicit none
        iadd = aux_loc + ivar- 1 + num_aux * (i - 1) + num_aux * mitot * (j - 1)
    end function iaddaux

    ! Output timing
    call system_clock(clock_start,clock_rate)
    call cpu_time(cpu_start)

    ! Count how many aux components requested
    output_aux_num = 0
    do i=1, num_aux
        output_aux_num = output_aux_num + output_aux_components(i)
    end do

    ! Note:  Currently outputs all aux components if any are requested
    out_aux = ((output_aux_num > 0) .and.               &
              ((.not. output_aux_onlyonce) .or. (time == t0)))
                
    ! Construct file names
    file_name(1) = 'fort.qxxxx'
    file_name(2) = 'fort.txxxx'
    file_name(3) = 'fort.axxxx'
    file_name(4) = 'fort.bxxxx'
    num_stop = frame
    do i = 10, 7, -1
        digit = mod(num_stop, 10)
        do j = 1, 4
            file_name(j)(i:i) = char(ichar('0') + digit)
        end do
        num_stop = num_stop / 10
    end do

    ! ==========================================================================
    ! Write out fort.q file (and fort.bXXXX or fort.bXXXX.nc files if necessary)
    ! Here we let fort.q be out_unit and the the other two be out_unit + 1
    open(unit=out_unit, file=file_name(1), status='unknown', form='formatted')
    if (output_style == 3) then
        open(unit=out_unit + 1, file=file_name(4), status="unknown",    &
             access='stream')
    else if (output_style == 4) then
        stop
    end if
    num_grids = 0

    ! Loop over levels
    do level = level_begin, level_end
        level_ptr = lstart(level)
        ! Loop over grids on each level
        do while (level_ptr /= 0)
            ! Extract grid data
            num_grids = num_grids + 1
            num_cells(1) = node(ndihi, level_ptr) - node(ndilo, level_ptr) + 1
            num_cells(2) = node(ndjhi, level_ptr) - node(ndjlo, level_ptr) + 1
            q_loc = node(store1, level_ptr)
            lower_corner = [rnode(cornxlo, level_ptr), rnode(cornylo, level_ptr)]
            delta = [hxposs(level), hyposs(level)]

            if (num_dim == 1) then
                write(out_unit, q_header_format_1d) level_ptr, level,          &
                                                    num_cells(1),              &
                                                    lower_corner(1),           &
                                                    delta(1)
            else
                write(out_unit, q_header_format_2d) level_ptr, level,          &
                                                    num_cells(1),              &
                                                    num_cells(2),              &
                                                    lower_corner(1),           &
                                                    lower_corner(2),           &
                                                    delta(1), delta(2)
            end if

            ! Output grids
            select case(output_format)
                ! ASCII output
                case(1)
                    ! Round off if nearly zero
                    do j = num_ghost + 1, num_cells(2) + num_ghost
                        do i = num_ghost + 1, num_cells(1) + num_ghost
                            write(out_unit, "(50e26.16)")                      &
                                            (alloc(iadd(m, i, j)), m=1, num_eqn)
                        end do
                        write(out_unit, *) ' '
                    end do

                ! What is case 2?
                case(2)
                    stop "Unknown format."

                ! Binary output
                case(3)
                    ! Note: We are writing out ghost cell data also
                    write(out_unit + 1) alloc(iadd(1, 1, 1):                   &
                                             (iadd(num_eqn,                    &
                                                 num_cells(1) + 2 * num_ghost, &
                                                 num_cells(2) + 2 * num_ghost))

                ! NetCDF output
                case(4)
#ifdef NETCDF
                    stop "NetCDF output not yet implemented!"
#else
                    print *, "ERROR:  NetCDF library is not available."
                    print *, "  Check the documentation as to how to include"
                    print *, "  the ability to output in NetCDF formats."
                    stop
#endif
                else
                    stop "Unsupported output format", output_format,"."

            end select
            level_ptr = node(levelptr, level_ptr)
        end do
    end do

    ! ==========================================================================
    ! Write out fort.a file
    stop "Aux files not implemented"

    ! ==========================================================================
    ! Write fort.t file
    open(unit=out_unit, file=file_name(2), status='unknown', form='formatted')

    ! Handle special case of using 2D AMR to do 1D AMR
    if (num_cells(2) > 1) then
        num_dim = 2
    else
        num_dum = 1
    end if

    ! Note:  We need to print out num_ghost too in order to strip ghost cells
    !        from q array when reading in pyclaw.io.binary
    write(out_unit, t_file_format) time, num_eqn, num_grids, num_aux, num_dim, &
                                   num_ghost
    close(out_unit)

    ! ==========================================================================
    ! Print output info
    print console_format, frame, time

    ! Increment frame counter
    frame = frame + 1

    ! Ouptut timing
    call system_clock(clock_finish,clock_rate)
    call cpu_time(cpu_finish)
    timeValout = timeValout + clock_finish - clock_start
    timeValoutCPU = timeValoutCPU + cpu_finish - cpu_start

end subroutine new_valout


