subroutine setprob

    implicit none

    real(kind=8) :: gamma
    common /cparam/  gamma

    ! ratio of specific heats, needed in Riemann solver:
    gamma = 1.4d0

end subroutine setprob
