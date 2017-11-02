module sweep_misc_module

    implicit none

    contains

attributes(device) &
subroutine compute_cqxx(cqxx, wave, s1, s2, dtdx)

    implicit none

    real(kind=8), intent(inout) :: cqxx(NEQNS)
    real(kind=8), intent(in) :: wave(NEQNS, NWAVES)
    real(kind=8), intent(in) :: s1, s2
    real(kind=8), intent(in) :: dtdx

    cqxx = 0.d0

    cqxx(1) = cqxx(1) + dabs(s1) * &
    (1.d0 - dabs(s1)*dtdx) * wave(1,1) &
    + dabs(s2) * &
    (1.d0 - dabs(s2)*dtdx) * wave(1,2)

    cqxx(2) = cqxx(2) +  dabs(s1) * &
    (1.d0 - dabs(s1)*dtdx) * wave(2,1) &
    +  dabs(s2) * &
    (1.d0 - dabs(s2)*dtdx) * wave(2,2)

    cqxx(3) = cqxx(3) +  dabs(s1) * &
    (1.d0 - dabs(s1)*dtdx) * wave(3,1) &
    +  dabs(s2) * &
    (1.d0 - dabs(s2)*dtdx) * wave(3,2)
end subroutine compute_cqxx

attributes(device) &
subroutine compute_cqyy(cqyy, wave, s1, s2, dtdy)

    implicit none

    real(kind=8), intent(inout) :: cqyy(NEQNS)
    real(kind=8), intent(in) :: wave(NEQNS, NWAVES)
    real(kind=8), intent(in) :: s1, s2
    real(kind=8), intent(in) :: dtdy

    cqyy = 0.d0

    cqyy(1) = cqyy(1) + dabs(s1) * &
    (1.d0 - dabs(s1)*dtdy) * wave(1,1) &
    + dabs(s2) * &
    (1.d0 - dabs(s2)*dtdy) * wave(1,2)

    cqyy(2) = cqyy(2) +  dabs(s1) * &
    (1.d0 - dabs(s1)*dtdy) * wave(2,1) &
    +  dabs(s2) * &
    (1.d0 - dabs(s2)*dtdy) * wave(2,2)

    cqyy(3) = cqyy(3) +  dabs(s1) * &
    (1.d0 - dabs(s1)*dtdy) * wave(3,1) &
    +  dabs(s2) * &
    (1.d0 - dabs(s2)*dtdy) * wave(3,2)
end subroutine compute_cqyy

end module sweep_misc_module
