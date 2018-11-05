module setprob_module

    save
    real(CLAW_REAL), parameter :: grav = 9.81, &
        dry_tolerance = 0.001, &
        rho = 1000.d0
#ifdef CUDA
    real(CLAW_REAL), constant :: grav_d, dry_tolerance_d, rho_d ! in device constant memory
#endif

    contains
    subroutine setprob
        implicit none

#ifdef CUDA
        grav_d = grav
        dry_tolerance_d = dry_tolerance
        rho_d = rho
#endif

        return
    end subroutine
end module
