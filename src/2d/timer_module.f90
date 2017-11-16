! =========================================================================
!  Subroutines for timing the code
!  1. cpu_timers is for timing wall time of a section of the code. To use 
!     it:
!     1) call take_cpu_timer(timer_name, timer_id). This will set timer with 
!        id **timer_id** to have name **timer_name**. This subroutine only 
!        takes effect on timer with id **timer_id** the first time it is called.
!     2) call cpu_timer_start(timer_id) with the same id **timer_id** at the 
!        beginning of the code section you want to time.
!     3) call cpu_timer_stop(timer_id) with the same id **timer_id** at the 
!        end of the code section you want to time.
!     4) A summary of all timers that have been taken by doing step 1) will be
!        printed out at the end of the program.
! =========================================================================
module timer_module
    implicit none
    save

    ! initialized to be 0
    integer, parameter :: max_cpu_timers = 20
    integer(kind=8) :: clock_rate

    type timer_type
        character(len=100) :: timer_name
        logical :: used = .false.
        logical :: running = .false.
        integer(kind=8) :: start_time = 0
        integer(kind=8) :: stop_time = 0
        integer(kind=8) :: accumulated_time = 0 ! in clock_cycle, not seconds
    end type timer_type


    ! measure wall time of CPU codes
    type(timer_type) :: cpu_timers(max_cpu_timers)
contains
    subroutine take_cpu_timer(timer_name_, timer_id)
        implicit none
        character(len=*), intent(in) :: timer_name_
        integer, intent(in) :: timer_id
        if (timer_id > max_cpu_timers) then
            print *, "timer_id for the cpu timer should be between 1 and ", max_cpu_timers
            stop
        endif
        if (.not. cpu_timers(timer_id)%used) then
            cpu_timers(timer_id)%timer_name = timer_name_
            cpu_timers(timer_id)%used = .true.
        else
            ! check to make sure we are not trying to give a new name to a existing timer
            if (cpu_timers(timer_id)%timer_name /= timer_name_) then
                print *, "Warning: trying to take a timer that's already assigned"
            endif
        endif
    end subroutine take_cpu_timer 

    subroutine cpu_timer_start(timer_id)
        implicit none
        integer, intent(in) :: timer_id
        if (.not. cpu_timers(timer_id)%used) then
            print *, "Warning: Trying to use a non-initialized cpu timer."
        endif
        if (.not. cpu_timers(timer_id)%running) then
            call system_clock(cpu_timers(timer_id)%start_time, clock_rate)
            cpu_timers(timer_id)%running = .true.
        else
            print *, "Warning: Trying to start a timer that's already running"
        endif
    end subroutine cpu_timer_start

    subroutine cpu_timer_stop(timer_id)
        implicit none
        integer, intent(in) :: timer_id
        if (.not. cpu_timers(timer_id)%used) then
            print *, "Warning: Trying to use a non-initialized cpu timer."
        endif
        if (cpu_timers(timer_id)%running) then
            call system_clock(cpu_timers(timer_id)%stop_time, clock_rate)
            cpu_timers(timer_id)%accumulated_time = cpu_timers(timer_id)%accumulated_time + &
                cpu_timers(timer_id)%stop_time - cpu_timers(timer_id)%start_time
            cpu_timers(timer_id)%running = .false.
        else
            print *, "Warning: Trying to stop a timer that's not running"
        endif
    end subroutine cpu_timer_stop

    subroutine print_all_cpu_timers()
        implicit none
        integer :: i
        print *, "Elapsed wall time recorded by all cpu timers: "
        do i = 1, max_cpu_timers
            if (cpu_timers(i)%used) then
                print *, "Wall time on ", cpu_timers(i)%timer_name, " is: ", &
                    real(cpu_timers(i)%accumulated_time,kind=8)/real(clock_rate,kind=8), " seconds."
            endif
        enddo
    end subroutine print_all_cpu_timers


end module timer_module
