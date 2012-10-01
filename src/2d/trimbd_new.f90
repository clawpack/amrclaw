!
! :::::::::::::::::::::::: TRIMBD :::::::::::::::::::::::::::;
!  if used array is completely set (=1.) then return set=true, 
!  otherwise return false, along with the dimensions of the smallest 
!  rectangle containing all unset points in il,ir,jb,jt.
! ::::::::::::::::::::::::::::::::::::;::::::::::::::::::::::;
!
subroutine trimbd(used,nrow,ncol,set,il,ir,jb,jt)

    implicit none

    ! Input
    integer, intent(in) :: nrow, ncol
    integer(kind=1), intent(in) :: used(nrow,ncol)

    ! Output
    logical, intent(out) :: set
    integer, intent(out) :: il, ir, jb, jt

    ! Locals
    integer :: i, j
    integer(kind=1) :: check

    if (sum(used) >= nrow * ncol ) then
        set = .true.
    else
        set = .false.
 
        check = 1
        do i = 1,nrow
             do j = 1,ncol
                check = min(check,used(i,j))
            enddo
            il = i
            if (check == 0) exit
        enddo

        check = 1
        do i = 1,nrow
            do j = 1,ncol
                check = min(check,used(nrow - i + 1,j))
            enddo
            ir = nrow - i + 1
            if (check == 0) exit
        enddo

        check = 1
        do j = 1,ncol
            do i = 1,nrow
                check = min(check,used(i,j))
            enddo
           jb = j
           if (check == 0) exit
        enddo
     
       check = 1
        do j = 1,ncol
            do i = 1,nrow
                check = min(check,used(i,ncol - j + 1))
            enddo
            jt = ncol - j + 1
            if (check == 0) exit
        enddo

    endif

end subroutine trimbd