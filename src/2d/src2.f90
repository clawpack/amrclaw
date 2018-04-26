subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)

!> Called to update q by solving source term equation 
!! $q_t = \psi(q)$ over time dt starting at time t.
!!
!! This default version does nothing. 
 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(CLAW_REAL), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(CLAW_REAL), intent(in) ::  aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(CLAW_REAL), intent(inout) ::  q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

end subroutine src2
