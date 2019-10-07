module setprob_module

    save
    real(CLAW_REAL) :: rho_l, rho_r, bulk_l, bulk_r


    contains
    subroutine setprob
        implicit real(CLAW_REAL) (a-h,o-z)
        character*25 fname

        !
        !     # Set the material parameters for the acoustic equations
        !
        !
        iunit = 7
        fname = 'setprob.data'
        !     # open the unit with new routine from Clawpack 4.4 to skip over
        !     # comment lines starting with #:
        call opendatafile(iunit, fname)

        !
        !     # Density and bulk modulus:

        read(7,*) rho_l
        read(7,*) rho_r
        read(7,*) bulk_l
        read(7,*) bulk_r

        return
    end subroutine
end module
