real(kind=8) function stream(x,y)
    ! ================================================
    !  Stream function in physical space (x,y).
    !  Clockwise rotation, rotates fully in time 1.

    implicit none

    real(kind=8), intent(in) :: x,y

    stream = 3.1415926535897931d0 *(x**2 + y**2)

end function stream
