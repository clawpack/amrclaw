!
! :::::::::::::::::::::::: TRIMBD :::::::::::::::::::::::::::;
!  if used array is completely set (=1.) then return set=true, 
!  otherwise return false, along with the dimensions of the smallest 
!  rectangle containing all unset points in il,ir.
! ::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::;
!
subroutine trimbd(used,nrow,set,unset_rect)

    implicit none

    ! Input
    integer, intent(in) :: nrow
    integer(kind=1), intent(in) :: used(nrow)

    ! Output
    logical, intent(out) :: set
    integer, intent(out) :: unset_rect(2)

    ! Locals
    integer :: i, utot
    integer(kind=1) :: check

       utot = 0
        do 100 i = 1,nrow
100        utot = utot + used(i)

    if (utot .eq. nrow ) then
        set = .true.
    else
        set = .false.
 
        check = 1
        do i = 1,nrow
            check = min(check,used(i))
            unset_rect(1) = i
            if (check == 0) exit
        enddo

        check = 1
        do i = 1,nrow
            check = min(check,used(nrow - i + 1))
            unset_rect(2) = nrow - i + 1
            if (check == 0) exit
        enddo

    endif

end subroutine trimbd
