subroutine setprob()

    use adjoint_module, only: read_adjoint_data

    implicit none

    call read_adjoint_data()    !# Read adjoint solution

end subroutine setprob
