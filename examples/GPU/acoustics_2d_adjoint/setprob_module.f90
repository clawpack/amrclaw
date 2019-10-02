module setprob_module

    save

    real(CLAW_REAL) :: rho, bulk, cc, zz
#ifdef CUDA
    real(CLAW_REAL), constant :: rho_d, bulk_d ! in device constant memory
#endif

    contains
subroutine setprob

    use adjoint_module, only: read_adjoint_data
    implicit none

    character*25 :: fname
    integer :: iunit
    
    real(CLAW_REAL) :: rho,bulk,cc,zz

    ! Set the material parameters for the acoustic equations
 
    iunit = 7
    fname = 'setprob.data'
    call opendatafile(iunit, fname)
              
    ! Density and bulk modulus:

    read(7,*) rho
    read(7,*) bulk

    ! Compute sound speed and impendance:

    cc = dsqrt(bulk/rho)
    zz = rho*cc

#ifdef CUDA
        rho_d = rho
        bulk_d = bulk
#endif

end subroutine setprob
end module
