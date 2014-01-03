real(kind=8) function psi(x,y)
    ! stream function 

    implicit none

    real(kind=8), intent(in) :: x, y

    real(kind=8) :: pi
    common /compsi/ pi

    psi = ((sin(pi * x))**2 * (sin(pi * y))**2) / pi

end function psi

