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

end module sweep_misc_module
