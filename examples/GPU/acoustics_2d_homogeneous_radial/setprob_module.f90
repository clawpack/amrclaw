module setprob_module

    save

    real(CLAW_REAL) :: rho, bulk, cc, zz
#ifdef CUDA
    real(CLAW_REAL), constant :: rho_d, bulk_d ! in device constant memory
#endif

    contains
    subroutine setprob
        implicit real(CLAW_REAL) (a-h,o-z)
        character*25 fname
        !
        !     # Set the material parameters for the acoustic equations
        !
        iunit = 7
        fname = 'setprob.data'
        !     # open the unit with new routine from Clawpack 4.4 to skip over
        !     # comment lines starting with #:
        call opendatafile(iunit, fname)
        !     # Density and bulk modulus:
        read(7,*) rho
        read(7,*) bulk
        !     # Compute sound speed and impendance:
        cc = sqrt(bulk/rho)
        zz = rho*cc

#ifdef CUDA
        rho_d = rho
        bulk_d = bulk
#endif

        return
    end subroutine
end module
