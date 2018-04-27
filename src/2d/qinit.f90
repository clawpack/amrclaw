! Default qinit file
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    !> Set initial conditions for the q array.
    ! This default version prints an error message since it should
    ! not be used directly.  Copy this to an application directory and
    ! loop over all grid cells to set values of q(1:meqn, 1:mx, 1:my).

    implicit none
    
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(CLAW_REAL), intent(in) :: xlower,ylower,dx,dy
    real(CLAW_REAL), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(CLAW_REAL), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    write(6,*) '*** Error -- you must set initial conditions'
    stop

end subroutine qinit
