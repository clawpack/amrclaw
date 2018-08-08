real(kind=8) pure function qtrue(x, y, t)

    implicit none

    ! Input
    real(kind=8), intent(in) :: x, y, t

    ! Locals
    real(kind=8) :: x0, y0, r

    ! Common block
    real(kind=8) :: ubar, vbar
    common /cparam/ ubar,vbar
      
    x0 = x - ubar * t
    y0 = y - vbar * t

    ! evaluate desired initial data at (x0,y0):
    r = sqrt((x0 + 0.2d0)**2 + (y0 - 0.4d0)**2)
    if (r <= 0.3d0) then
        qtrue = 1.d0
    else
        qtrue = 0.d0
    endif

    r = sqrt((x0 - 0.3d0)**2 + (y0 - 0.1d0)**2)
    qtrue = qtrue - exp(-15.d0 * r**2)

end function qtrue