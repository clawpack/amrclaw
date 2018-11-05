!> Output the results for a general system of conservation laws
!! in 2 dimensions
!!
!! Write the results to the file fort.q<iframe>
!! Use format required by matlab script  plotclaw2.m or Python tools
!!
!! set outaux = .true. to also output the aux arrays to fort.a<iframe>
subroutine valout(level_begin, level_end, time, num_eqn, num_aux)

    use amr_module, only: alloc, t0, output_aux_onlyonce, output_aux_components
    use amr_module, only: frame => matlabu, num_ghost => nghost, lstart
    use amr_module, only: hxposs, output_format, store1, storeaux
    use amr_module, only: node, rnode, ndilo, ndihi, cornxlo, mxnest, levelptr
    use amr_module, only: timeValout, timeValoutCPU, tvoll, tvollCPU, rvoll
    use amr_module, only: timeTick, tick_clock_start, t0

#ifdef HDF5
    use hdf5
#endif

    implicit none

    ! Input
    integer, intent(in) :: level_begin, level_end, num_eqn, num_aux
    real(kind=8), intent(in) :: time

    ! Locals
    logical :: timing_file_exists
    integer, parameter :: out_unit = 50
    integer :: i, j, m, level, output_aux_num, num_stop, digit
    integer :: grid_ptr, num_cells, num_grids, q_loc, aux_loc
    real(kind=8) :: lower_corner, delta
    logical :: out_aux
    character(len=10) :: file_name(5)

    ! Timing
    integer :: clock_start, clock_finish, clock_rate
    integer    tick_clock_finish, tick_clock_rate, timeTick_int
    real(kind=8) :: cpu_start, cpu_finish, t_CPU_overall, timeTick_overall
    character(len=256) :: timing_line, timing_substr
    character(len=*), parameter :: timing_file_name = "timing.csv"

    character(len=*), parameter :: header_format =                             &
                                    "(i6,'                 grid_number',/," // &
                                     "i6,'                 AMR_level',/,"   // &
                                     "i6,'                 mx',/,"          // &
                                     "e26.16,'    xlow', /, "               // &
                                     "e26.16,'    dx', /)"
    character(len=*), parameter :: t_file_format = "(e18.8,'    time', /,"  // &
                                           "i6,'                 meqn'/,"   // &
                                           "i6,'                 ngrids'/," // &
                                           "i6,'                 naux'/,"   // &
                                           "i6,'                 ndim'/,"   // &
                                           "i6,'                 nghost'/,/)"
    character(len=*), parameter :: timing_header_format =                      &
                                                  "(' wall time (', i2,')," // &
                                                  " CPU time (', i2,'), "   // &
                                                  "cells updated (', i2,'),')"
    character(len=*), parameter :: console_format = &
             "('AMRCLAW: Frame ',i4,' output files done at time t = ', d13.6,/)"

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
              ((.not. output_aux_onlyonce) .or. (abs(time - t0) < 1d-90)))
                
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
    if (output_format == 3) then
        open(unit=out_unit + 1, file=file_name(4), status="unknown",    &
             access='stream')
    else if (output_format == 4) then
        stop
        ! open(unit=out_unit + 1, file=file_name(4) // ".nc", status="unknown",  &
        !      access='stream')
    end if
    num_grids = 0

    ! Loop over levels
    do level = level_begin, level_end
        grid_ptr = lstart(level)
        delta = hxposs(level)

        ! Loop over grids on each level
        do while (grid_ptr /= 0)
            ! Extract grid data
            num_grids = num_grids + 1
            num_cells = node(ndihi, grid_ptr) - node(ndilo, grid_ptr) + 1
            q_loc = node(store1, grid_ptr)
            lower_corner = rnode(cornxlo, grid_ptr)

            ! Write out header data                 
            write(out_unit, header_format) grid_ptr, level, num_cells,         &
                                           lower_corner, delta
            ! Output grids
            select case(output_format)
                ! ASCII output
                case(1)
                    ! Round off if nearly zero
                    forall (m = 1:num_eqn,                              &
                            i=num_ghost + 1:num_cells + num_ghost,   &
                            abs(alloc(iadd(m, i))) < 1d-90)

                        alloc(iadd(m, i)) = 0.d0
                    end forall

                    do i = num_ghost + 1, num_cells + num_ghost
                        write(out_unit, "(50e26.16)")                      &
                                            (alloc(iadd(m, i)), m=1, num_eqn)
                    end do
                    write(out_unit, *) ' '

                ! What is case 2?
                case(2)
                    stop "Unknown format."

                ! Binary output
                case(3)
                    ! Note: We are writing out ghost cell data also
                    i = iadd(num_eqn, num_cells + 2 * num_ghost)
                    write(out_unit + 1) alloc(iadd(1, 1):i)

                ! HDF5 output
                case(4)
#ifdef HDF5
                    stop "HDF5 output not yet implemented!"
#else
                    print *, "ERROR:  HDF5 library is not available."
                    print *, "  Check the documentation as to how to include"
                    print *, "  the ability to output in HDF5 formats."
                    stop
#endif
                case default
                    print *, "Unsupported output format", output_format,"."
                    stop 

            end select
            grid_ptr = node(levelptr, grid_ptr)
        end do
    end do

    ! ==========================================================================
    ! Write out fort.a file
    if (out_aux) then
        if (output_format == 1) then
            open(unit=out_unit, file=file_name(3), status='unknown',        &
                 form='formatted')
        else if (output_format == 3) then
            open(unit=out_unit, file=file_name(3), status='unknown',        &
                 access='stream')
        else if (output_format == 4) then
            stop
            ! open(unit=out_unit, file=file_name(3) // ".nc", status='unknown',  &
            !      access='stream')
        end if

        do level = level_begin, level_end
            grid_ptr = lstart(level)
            delta = hxposs(level)

            ! Loop over grids on each level
            do while (grid_ptr /= 0)
                ! Extract grid data
                num_cells = node(ndihi, grid_ptr) - node(ndilo, grid_ptr) + 1
                aux_loc = node(storeaux, grid_ptr)
                lower_corner = rnode(cornxlo, grid_ptr)

                ! Output grids
                select case(output_format)
                    ! ASCII output
                    case(1)

                        ! We only output header info for aux data if writing 
                        ! ASCII data
                        write(out_unit, header_format) grid_ptr, level,     &
                                                       num_cells,           &
                                                       lower_corner,        &
                                                       delta

                        ! Round off if nearly zero
                        forall (m = 1:num_aux,                              &
                                i=num_ghost + 1:num_cells + num_ghost,      &
                                abs(alloc(iaddaux(m, i))) < 1d-90)

                            alloc(iaddaux(m, i)) = 0.d0
                        end forall

                        do i = num_ghost + 1, num_cells + num_ghost
                            write(out_unit, "(50e26.16)")                   &
                                         (alloc(iaddaux(m, i)), m=1, num_aux)
                        end do
                        write(out_unit, *) ' '

                    ! What is case 2?
                    case(2)
                        stop "Unknown format."

                    ! Binary output
                    case(3)
                        ! Note: We are writing out ghost cell data also
                        i = iaddaux(num_aux, num_cells + 2 * num_ghost)
                        write(out_unit) alloc(iaddaux(1, 1):i)

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
                    case default
                        print *, "Unsupported output format", output_format,"."
                        stop 

                end select
                grid_ptr = node(levelptr, grid_ptr)
            end do
        end do
    end if

    ! ==========================================================================
    ! Write fort.t file
    open(unit=out_unit, file=file_name(2), status='unknown', form='formatted')

    ! Note:  We need to print out num_ghost too in order to strip ghost cells
    !        from q array when reading in pyclaw.io.binary
    write(out_unit, t_file_format) time, num_eqn, num_grids, num_aux, 1,    &
                                   num_ghost
    close(out_unit)

    ! ==========================================================================
    ! Write out timing stats
    ! Assume that this has been started some where
    open(unit=out_unit, file=timing_file_name, form='formatted',         &
         status='old', action='write', position='append')
    
    timing_line = "(e16.6, ', ', e16.6, ', ', e16.6,"
    do level=1, mxnest
        timing_substr = "', ', e16.6, ', ', e16.6, ', ', e16.6"
        timing_line = trim(timing_line) // timing_substr
    end do
    timing_line = trim(timing_line) // ")"

    if (time == t0) then
        t_CPU_overall = 0.d0
        timeTick_overall = 0.d0
      else
        call cpu_time(t_CPU_overall)
        call system_clock(tick_clock_finish,tick_clock_rate)
        ! not including timeTick since not yet set up for restarting
        timeTick_int = tick_clock_finish - tick_clock_start
        timeTick_overall = real(timeTick_int, kind=8)/real(tick_clock_rate,kind=8)
      endif

    write(out_unit, timing_line) time, timeTick_overall, t_CPU_overall, &
        (real(tvoll(i), kind=8) / real(clock_rate, kind=8), &
         tvollCPU(i), rvoll(i), i=1,mxnest)
    
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

contains

    ! Index into q array
    pure integer function iadd(m, i)
        implicit none
        integer, intent(in) :: m, i
        iadd = q_loc + m - 1 + num_eqn * (i - 1)
    end function iadd

    ! Index into aux array
    pure integer function iaddaux(m, i)
        implicit none
        integer, intent(in) :: m, i
        iaddaux = aux_loc + m - 1 + num_aux * (i - 1)
    end function iaddaux

end subroutine valout


