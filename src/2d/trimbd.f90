!
! :::::::::::::::::::::::: TRIMBD :::::::::::::::::::::::::::;
!  if used array is completely set (=1.) then return set=true, 
!  otherwise return false, along with the dimensions of the smallest 
!  rectangle containing all unset points in il,ir,jb,jt.
! ::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::;
!
subroutine trimbd(used,nrow,ncol,set,unset_rect)

    implicit none

    ! Input
    integer, intent(in) :: nrow, ncol
    integer(kind=1), intent(in) :: used(nrow,ncol)

    ! Output
    logical, intent(out) :: set
!     integer, intent(out) :: il, ir, jb, jt
    integer, intent(out) :: unset_rect(4)

    ! Locals
    integer :: i, j, utot
    integer(kind=1) :: check

       utot = 0
        do 100 j = 1,ncol
        do 100 i = 1,nrow
100        utot = utot + used(i,j)

    if (utot .eq. nrow * ncol ) then
        set = .true.
    else
        set = .false.
 
        check = 1
        do i = 1,nrow
             do j = 1,ncol
                check = min(check,used(i,j))
            enddo
            unset_rect(1) = i
            if (check == 0) exit
        enddo

        check = 1
        do i = 1,nrow
            do j = 1,ncol
                check = min(check,used(nrow - i + 1,j))
            enddo
            unset_rect(2) = nrow - i + 1
            if (check == 0) exit
        enddo

        check = 1
        do j = 1,ncol
            do i = 1,nrow
                check = min(check,used(i,j))
            enddo
           unset_rect(3) = j
           if (check == 0) exit
        enddo
     
       check = 1
        do j = 1,ncol
            do i = 1,nrow
                check = min(check,used(i,ncol - j + 1))
            enddo
            unset_rect(4) = ncol - j + 1
            if (check == 0) exit
        enddo

    endif

end subroutine trimbd
