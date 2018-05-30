subroutine setprob

    use adjoint_module, only: read_adjoint_data
    implicit none

    character*25 :: fname
    integer :: iunit
    
    real(kind=8) :: rho,bulk,cc,zz
    common /cparam/ rho,bulk,cc,zz

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

end subroutine setprob
