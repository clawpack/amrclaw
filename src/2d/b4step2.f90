!> Called before each call to step2.
!! Use to set time-dependent aux arrays or perform other tasks.
!!
!! This default version does nothing. 
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(CLAW_REAL), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(CLAW_REAL), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(CLAW_REAL), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

end subroutine b4step2
